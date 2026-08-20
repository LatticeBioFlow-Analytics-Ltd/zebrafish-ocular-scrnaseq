"""Figure generation and end-to-end pipeline orchestration."""

from __future__ import annotations

import logging
import textwrap
from contextlib import contextmanager
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from . import figures as figs
from .cluster import process, subcluster
from .compartments import EVIDENCE_COLUMNS, assign_compartments, subset_compartment
from .config import Config
from .provenance import write_manifest
from .qc import build_qc_dataset

log = logging.getLogger(__name__)

# Inter-panel spacing for scanpy's own multi-panel grids (feature plots with a
# colourbar per panel). Categorical panels no longer rely on this: their
# legends are laid out explicitly in the embedding helpers below, because a
# spacing value wide enough for a 30-entry legend wastes half the page on every
# figure that lacks one.
PANEL_SPACING = 0.30


@contextmanager
def headless_plotting():
    """Render figures to file without a display, restoring the caller's backend.

    The backend is deliberately NOT set at import time. `screye/__init__.py`
    imports this module, so an import-time `matplotlib.use("Agg")` would also
    apply to the notebook: scanpy's `savefig_or_show` defaults `show=None` to
    `settings.autoshow`, calls `plt.show()`, and matplotlib warns
    "FigureCanvasAgg is non-interactive" for every inline figure.

    Under an inline backend the switch is skipped entirely, so calling `run()`
    from a notebook leaves interactive plotting intact.
    """
    previous = matplotlib.get_backend()
    interactive = "inline" in previous.lower() or "ipympl" in previous.lower()
    if not interactive:
        plt.switch_backend("Agg")
    try:
        yield
    finally:
        if not interactive and matplotlib.get_backend() != previous:
            plt.switch_backend(previous)


def _save(fig_dir: Path, name: str, width: float | None = figs.DOUBLE_COL) -> None:
    """Save a figure at journal width. Thin wrapper so every figure in the
    project goes through one place and inherits the same size, dpi and formats."""
    figs.save_figure(fig_dir, name, width=width)


# Dot sizing. scanpy's default largest dot is 200 pt^2; at the figs.DOT_PITCH
# grid (3.3 mm) a dot may be at most ~9 pt across before it touches its
# neighbour, so the ceiling is sqrt(45) ~ 6.7 pt.
LARGEST_DOT = 45.0
SMALLEST_DOT = 1.0
DOT_EDGE_LW = 0.25


def _sample_timepoints(adata: ad.AnnData) -> dict[str, str]:
    """Ordered mapping of sample name -> timepoint label."""
    obs = adata.obs
    return {
        str(s): str(obs.loc[obs["sample"] == s, "timepoint"].iloc[0])
        for s in obs["sample"].cat.categories
        if (obs["sample"] == s).any()
    }


# --- dotplots ----------------------------------------------------------------

def _chunk_var_groups(groups: dict[str, list[str]],
                      max_cols: int = figs.MAX_DOT_COLUMNS) -> list[dict[str, list[str]]]:
    """Split a var-group mapping into chunks that fit a double column.

    A dotplot column costs ~3.3 mm at 8 pt (see figs.DOT_PITCH), so only
    ~MAX_DOT_COLUMNS genes fit a 180 mm figure. Anything wider must become
    several figures: shrinking instead leaves the fixed-point dots overlapping
    and the rotated gene labels on top of each other. Groups are kept intact
    within a chunk; a single group wider than the budget is split and its
    parts numbered.
    """
    flat: dict[str, list[str]] = {}
    for name, genes in groups.items():
        if len(genes) <= max_cols:
            flat[name] = genes
        else:
            n_parts = -(-len(genes) // max_cols)
            per = -(-len(genes) // n_parts)
            for i in range(n_parts):
                flat[f"{name} ({i + 1}/{n_parts})"] = genes[i * per:(i + 1) * per]

    chunks: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] = {}
    count = 0
    for name, genes in flat.items():
        if current and count + len(genes) > max_cols:
            chunks.append(current)
            current, count = {}, 0
        current[name] = genes
        count += len(genes)
    if current:
        chunks.append(current)
    return chunks


def _style_dotplot(plot, fig_dir: Path, name: str) -> None:
    """Apply shared pitch, dot sizing and padding, then save at natural size.

    The per-category pitch is set to figs.DOT_PITCH (scanpy's default of
    0.37"/column puts a 60-gene panel at 565 mm) and scanpy's own layout is
    then left to size the figure from its rows and columns. Saved with
    `width=None`: rescaling a grid of fixed-point dots and labels to a column
    width is what produced the overlaps; the column budget is enforced upstream
    by _chunk_var_groups instead. save_figure still measures the rendered width
    and warns if it exceeds the page.
    """
    plot.DEFAULT_CATEGORY_WIDTH = figs.DOT_PITCH
    plot.DEFAULT_CATEGORY_HEIGHT = figs.DOT_PITCH
    plot.style(
        largest_dot=LARGEST_DOT,
        smallest_dot=SMALLEST_DOT,
        dot_edge_lw=DOT_EDGE_LW,
        x_padding=0.45,
        y_padding=0.45,
    ).legend(width=figs.DOT_LEGEND_WIDTH)  # default 1.5 overlaps the size key
    plot.make_figure()
    figs.save_figure(fig_dir, name, width=None)


def _dotplot_chunked(adata: ad.AnnData, groups: dict[str, list[str]],
                     groupby: str, fig_dir: Path, name: str,
                     dendrogram: bool = False) -> None:
    """Draw one content-sized dotplot per width chunk.

    Every chunk keeps all `groupby` categories as rows, so cross-cluster
    comparison survives the split - only the gene columns are divided.
    """
    chunks = _chunk_var_groups(groups)
    if len(chunks) > 1:
        log.info("%s: %d gene columns split across %d figures of <=%d columns",
                 name, sum(len(v) for v in groups.values()), len(chunks),
                 figs.MAX_DOT_COLUMNS)
    # The chunk count depends on the data (cluster count, genes per set), so a
    # re-run can produce fewer files than the last one - or one where the last
    # produced several. Remove every output under this name first, or stale
    # figures from the previous run sit beside fresh ones indistinguishably.
    for stale in list(fig_dir.glob(f"{name}.*")) + list(fig_dir.glob(f"{name}_[0-9]*.*")):
        stale.unlink()
    for i, chunk in enumerate(chunks, start=1):
        plot = sc.pl.dotplot(
            adata, chunk, groupby=groupby, use_raw=True, standard_scale="var",
            dendrogram=dendrogram, return_fig=True, show=False,
        )
        suffix = "" if len(chunks) == 1 else f"_{i}"
        _style_dotplot(plot, fig_dir, f"{name}{suffix}")


def _de_dotplot(adata: ad.AnnData, fig_dir: Path, name: str,
                n_genes: int = 4, groupby: str = "leiden") -> None:
    """Top differentially expressed genes per cluster, as sized dotplots."""
    ranked = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])

    # Order the per-cluster gene blocks by the expression dendrogram so related
    # clusters sit adjacently and the diagonal block structure is readable.
    dendro_key = f"dendrogram_{groupby}"
    if dendro_key not in adata.uns:
        sc.tl.dendrogram(adata, groupby=groupby)
    order = [str(g) for g in adata.uns[dendro_key]["categories_ordered"]]

    top = {
        group: list(dict.fromkeys(ranked[group].iloc[:n_genes]))
        for group in order if group in ranked.columns
    }
    _dotplot_chunked(adata, top, groupby, fig_dir, name, dendrogram=True)


def _marker_dotplot(adata: ad.AnnData, markers, groupby: str,
                    fig_dir: Path, name: str) -> None:
    """Curated marker sets or genes of interest, as sized dotplots."""
    if isinstance(markers, dict):
        # Multi-line group labels: a 30-character set name rotated 90 degrees
        # is a 50 mm tower that dwarfs the data. Wrapped at ~14 characters it
        # is a compact block at the true 8 pt size.
        display = {
            textwrap.fill(k.replace("_", " "), width=14): v
            for k, v in markers.items()
        }
        _dotplot_chunked(adata, display, groupby, fig_dir, name)
        return

    genes = list(markers)
    plot = sc.pl.dotplot(
        adata, genes, groupby=groupby, use_raw=True, standard_scale="var",
        return_fig=True, show=False,
    )
    _style_dotplot(plot, fig_dir, name)


# --- embeddings ---------------------------------------------------------------

def _panel_title(key: str) -> str:
    return key.replace("_", " ")


def _palette_for(adata: ad.AnnData, key: str) -> list[str] | None:
    """Okabe-Ito for keys with few enough levels; scanpy's default otherwise."""
    n = adata.obs[key].nunique()
    return figs.CATEGORICAL if n <= len(figs.CATEGORICAL) else None


def _draw_embedding(adata: ad.AnnData, basis: str, key: str, ax,
                    on_data: bool = False) -> None:
    """One embedding panel on a caller-owned axis, legend handled by caller."""
    sc.pl.embedding(
        adata, basis=basis, color=key, ax=ax, show=False, frameon=False,
        title=_panel_title(key), palette=_palette_for(adata, key),
        legend_loc="on data" if on_data else "none",
        legend_fontsize=figs.LABEL_PT, legend_fontoutline=2,
    )


def _legend_handles(adata: ad.AnnData, key: str) -> list[mlines.Line2D]:
    cats = list(adata.obs[key].cat.categories)
    colours = list(adata.uns.get(f"{key}_colors", figs.CATEGORICAL))
    return [
        mlines.Line2D([], [], marker="o", linestyle="", markersize=4,
                      markeredgewidth=0, color=c, label=str(cat).replace("_", " "))
        for cat, c in zip(cats, colours)
    ]


def _embedding_overview(adata: ad.AnnData, basis: str, fig_dir: Path,
                        name: str) -> None:
    """Cluster and timepoint panels side by side.

    Clusters are labelled on the data rather than in a legend: a 30-cluster
    legend is four columns of swatches that no spacing value can accommodate at
    8 pt within a 180 mm figure, and matching 30 colours between a panel and a
    distant key is not something a reader can actually do. The timepoint legend
    has few entries and sits beneath its panel.
    """
    panel = (figs.DOUBLE_COL - 0.4) / 2
    _fig, axes = plt.subplots(1, 2, figsize=(figs.DOUBLE_COL, panel + 0.45))
    _draw_embedding(adata, basis, "leiden", axes[0], on_data=True)
    _draw_embedding(adata, basis, "timepoint", axes[1])
    axes[1].legend(
        handles=_legend_handles(adata, "timepoint"),
        loc="upper center", bbox_to_anchor=(0.5, -0.02),
        ncols=min(4, adata.obs["timepoint"].nunique()),
        frameon=False, handletextpad=0.3, columnspacing=1.0,
    )
    _save(fig_dir, name)


def _embedding_side_legend(adata: ad.AnnData, basis: str, key: str,
                           fig_dir: Path, name: str,
                           max_legend: int = 48) -> None:
    """One panel with its legend laid out in reserved space to the right.

    The panel takes ~55% of the width and the legend the rest, so a long key
    ("periocular_mesenchyme (low confidence)") has room instead of being drawn
    over the neighbouring panel. Beyond `max_legend` entries a key stops being
    readable as a legend at all; the panel is still drawn and the omission
    logged.
    """
    n = adata.obs[key].nunique()
    fig = plt.figure(figsize=(figs.DOUBLE_COL, 0.52 * figs.DOUBLE_COL))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.55, 0.45], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    fig.add_subplot(gs[0, 1]).set_axis_off()  # reserved for the legend
    _draw_embedding(adata, basis, key, ax)
    if n <= max_legend:
        ax.legend(
            handles=_legend_handles(adata, key),
            loc="center left", bbox_to_anchor=(1.02, 0.5),
            ncols=1 + (n - 1) // 24, frameon=False,
            handletextpad=0.3, labelspacing=0.4, borderaxespad=0,
        )
    else:
        log.info("%s: %d levels in %r is too many for a readable legend; "
                 "see the composition tables instead", name, n, key)
    _save(fig_dir, name)


def plot_embeddings(adata: ad.AnnData, fig_dir: Path, compute_tsne: bool) -> None:
    """UMAP and t-SNE coloured by cluster, timepoint, cell type and library."""
    bases = ["umap"]
    if compute_tsne and "X_tsne" in adata.obsm:
        bases.append("tsne")

    sample_to_tp = _sample_timepoints(adata)
    # One library per timepoint makes `sample` a recolouring of `timepoint`;
    # a separate panel would claim information that is not there.
    sample_redundant = len(sample_to_tp) == len(set(sample_to_tp.values()))

    for basis in bases:
        _embedding_overview(adata, basis, fig_dir, f"{basis}_overview")
        _embedding_side_legend(adata, basis, "cell_type", fig_dir,
                               f"{basis}_cell_types")
        if not sample_redundant:
            _embedding_side_legend(adata, basis, "sample", fig_dir,
                                   f"{basis}_by_library")

    # Sample composition per cluster: a cluster drawn from a single library is
    # either stage-specific biology or an uncorrected technical artefact, and
    # the two are worth telling apart before the labels are trusted. Libraries
    # are shaded by timepoint so both levels stay readable with many libraries.
    shades = figs.timepoint_shades(sample_to_tp)
    comp = (
        pd.crosstab(adata.obs["leiden"], adata.obs["sample"], normalize="index")
        .sort_index()
    )
    _fig, ax = figs.figure(width=figs.DOUBLE_COL, height=0.45 * figs.DOUBLE_COL)
    comp.plot(kind="bar", stacked=True, width=0.85, ax=ax, legend=False,
              color=[shades[str(s)] for s in comp.columns])
    ax.set_ylabel("fraction of cluster")
    ax.set_xlabel("leiden cluster")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title=None, bbox_to_anchor=(1.02, 1), loc="upper left",
              ncols=1 + (len(comp.columns) - 1) // 14, frameon=False)
    _save(fig_dir, "cluster_composition_by_sample")


def plot_markers(adata: ad.AnnData, markers: dict[str, list[str]], fig_dir: Path) -> None:
    """Dotplot of the curated panel and the data-driven top genes."""
    present = {
        name: [g for g in genes if g in adata.raw.var_names]
        for name, genes in markers.items()
    }
    present = {k: v for k, v in present.items() if v}
    if present:
        _marker_dotplot(adata, present, "leiden", fig_dir, "dotplot_curated_markers")

    _de_dotplot(adata, fig_dir, "dotplot_top_de_genes")


def plot_qc(adata: ad.AnnData, fig_dir: Path) -> None:
    """Per-sample QC distributions, one metric per row.

    Rows rather than scanpy's side-by-side panels so that the width available
    to each metric is the full figure and scales with the number of libraries:
    three panels of 22 violins each cannot share 180 mm legibly.
    """
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    n_samples = adata.obs["sample"].nunique()
    width = figs.SINGLE_COL if n_samples <= 4 else figs.DOUBLE_COL

    # Same colour logic as every other per-library figure: libraries shaded by
    # timepoint. Left to seaborn, day5 renders blue and day8 orange - the exact
    # reverse of the Okabe-Ito assignment used everywhere else in the project.
    shades = figs.timepoint_shades(_sample_timepoints(adata))
    palette = [shades[str(s)] for s in adata.obs["sample"].cat.categories]

    _fig, axes = plt.subplots(
        len(metrics), 1, figsize=(width, 1.35 * len(metrics)), sharex=True,
    )
    for ax, metric in zip(np.atleast_1d(axes), metrics):
        # dodge=False is required. With seaborn >= 0.13 scanpy passes
        # hue=groupby to sns.violinplot but not to the stripplot beneath, so a
        # dodged violin and its points no longer overlay. The strip itself is
        # off: at these cell counts it renders as a solid block that hides the
        # distribution. inner="box" substitutes an interquartile box, which is
        # also what DMM asks for when points are too numerous to plot.
        # hue="sample" (the groupby itself) with legend=False is seaborn 0.13's
        # required idiom for per-category palettes; palette without hue is
        # deprecated and warns as FutureWarning, which the test suite escalates.
        sc.pl.violin(
            adata, metric, groupby="sample", show=False, ax=ax,
            rotation=45 if n_samples <= 8 else 90, palette=palette,
            hue="sample", legend=False,
            dodge=False, stripplot=False, inner="box", density_norm="width",
        )
    for ax in np.atleast_1d(axes)[:-1]:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)
    _save(fig_dir, "qc_violin_by_sample", width=width)
    (fig_dir / "qc_violin_by_sample_legend.txt").write_text(
        "Per-library distributions of detected genes, total UMI counts and "
        "mitochondrial read percentage, before filtering. Violins show the full "
        "kernel density; the internal box marks the interquartile range with "
        "whiskers extending to 1.5x the IQR and the white point at the median.\n"
    )

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", color="sample", show=False)
    _save(fig_dir, "qc_counts_vs_mito")


def plot_genes_of_interest(
    adata: ad.AnnData,
    genes: tuple[str, ...] | list[str],
    fig_dir: Path,
    groupby: str = "cell_type",
) -> pd.DataFrame:
    """Plot named genes across the embedding, clusters and cell types.

    Returns a tidy table of mean expression and the fraction of cells with any
    detected transcript, per gene per group. The table matters as much as the
    figures: a gene can look convincingly present on a UMAP while being detected
    in 3% of the cells it appears to mark, and colour scales do not show that.

    Genes are read from `.raw`, which holds the full log-normalised matrix. The
    analysis object is subset to highly variable genes, and a transporter
    expressed narrowly in one cell type is often not highly variable across the
    whole dataset, so looking it up in `.var_names` would report it as absent.
    """
    if not genes:
        return pd.DataFrame()

    source = adata.raw
    present = [g for g in genes if g in source.var_names]
    missing = [g for g in genes if g not in source.var_names]
    if missing:
        log.warning(
            "Not found in the data: %s. Check the symbol against ZFIN and "
            "against the CellRanger reference GTF before concluding the gene "
            "is unexpressed.", ", ".join(missing),
        )
    if not present:
        log.warning("None of the requested genes are present; skipping.")
        return pd.DataFrame()

    # Separate "in the reference but never detected" from "detected somewhere".
    # A gene with no non-zero counts still plots: scanpy falls back to a
    # symmetric +/-0.1 colour scale, which looks like real, centred variation
    # and is one of the easier ways to read structure into an empty panel.
    detected, undetected = [], []
    for gene in present:
        column = source[:, gene].X
        total = column.sum() if not hasattr(column, "nnz") else column.nnz
        (detected if total > 0 else undetected).append(gene)

    if undetected:
        log.warning(
            "Present in the reference but not detected in any cell: %s. "
            "Excluded from the expression panels, since an all-zero gene is "
            "rendered on an artificial colour scale that resembles real "
            "variation. Reported in the summary table with zero counts.",
            ", ".join(undetected),
        )
    if not detected:
        log.warning("No requested gene is detected in any cell; skipping panels.")
        return pd.DataFrame()

    # 1. Where in the embedding. vmax="p99" stops a handful of high-count cells
    #    flattening the colour scale for everything else.
    figs.scanpy_panel_size(n_cols=3)
    sc.pl.umap(adata, color=detected, use_raw=True, ncols=3, cmap=figs.SEQUENTIAL,
               vmin=0, vmax="p99", frameon=False, wspace=PANEL_SPACING, show=False)
    _save(fig_dir, "goi_umap")

    if "X_tsne" in adata.obsm:
        figs.scanpy_panel_size(n_cols=3)
        sc.pl.tsne(adata, color=detected, use_raw=True, ncols=3, cmap=figs.SEQUENTIAL,
                   vmin=0, vmax="p99", frameon=False, wspace=PANEL_SPACING, show=False)
        _save(fig_dir, "goi_tsne")

    # 2. Dot size encodes the fraction of cells expressing, which is the part a
    #    feature plot cannot show.
    for key, tag in ((groupby, "celltype"), ("leiden", "cluster")):
        if key in adata.obs:
            _marker_dotplot(adata, detected, key, fig_dir, f"goi_dotplot_by_{tag}")

    # 3. Distribution per cell type, split by timepoint.
    sc.pl.stacked_violin(adata, detected, groupby=groupby, use_raw=True,
                         standard_scale="var", show=False)
    _save(fig_dir, "goi_stacked_violin_by_celltype")

    n_samples = adata.obs["sample"].nunique()
    for gene in detected:
        _fig, ax = figs.figure(
            width=figs.SINGLE_COL if n_samples <= 4 else figs.DOUBLE_COL,
            height=1.6,
        )
        sc.pl.violin(adata, gene, groupby="sample", use_raw=True, ax=ax,
                     rotation=45 if n_samples <= 8 else 90,
                     dodge=False, stripplot=False, inner="box",
                     density_norm="width", show=False)
        _save(fig_dir, f"goi_violin_{gene}_by_sample",
              width=figs.SINGLE_COL if n_samples <= 4 else figs.DOUBLE_COL)

    # 4. The numbers behind the pictures.
    matrix = source[:, present].X
    values = pd.DataFrame(
        matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix),
        index=adata.obs_names, columns=present,
    )
    values[[groupby, "sample"]] = adata.obs[[groupby, "sample"]].astype(str).to_numpy()

    grouped = values.groupby([groupby, "sample"], observed=True)
    summary = (
        grouped[present].mean().stack().rename("mean_expression").reset_index()
        .rename(columns={"level_2": "gene"})
        .merge(
            grouped[present].apply(lambda d: (d > 0).mean()).stack()
            .rename("fraction_expressing").reset_index()
            .rename(columns={"level_2": "gene"}),
            on=[groupby, "sample", "gene"],
        )
    )
    summary["n_cells"] = summary.apply(
        lambda r: int(((values[groupby] == r[groupby]) & (values["sample"] == r["sample"])).sum()),
        axis=1,
    )
    return summary.sort_values(["gene", groupby, "sample"]).reset_index(drop=True)


def plot_compartments(adata: ad.AnnData, evidence: pd.DataFrame, fig_dir: Path) -> None:
    """Pass-1 overview: where each compartment sits, and the evidence for it."""
    _embedding_side_legend(adata, "umap", "compartment_group", fig_dir,
                           "pass1_compartment_umap")
    _embedding_side_legend(adata, "umap", "compartment", fig_dir,
                           "pass1_compartment_fine_umap")

    # Score heatmap. Diverging red/blue centred on zero: a compartment score is
    # a difference from a matched random reference set, so zero is the
    # meaningful midpoint and a sequential scale would misrepresent it.
    score_cols = [c for c in evidence.columns if c not in EVIDENCE_COLUMNS]
    matrix = evidence[score_cols]
    # Row pitch as in the dotplots: enough for an 8 pt row label, so the figure
    # grows with the cluster count instead of squeezing rows together.
    height = min(max(1.6, figs.DOT_PITCH * len(matrix) + 1.1), figs.WORKING_HEIGHT)
    fig, ax = figs.figure(width=figs.DOUBLE_COL, height=height)
    vmax = float(np.nanmax(np.abs(matrix.to_numpy())))
    im = ax.imshow(matrix.to_numpy(), aspect="auto", cmap=figs.DIVERGING,
                   norm=figs.diverging_norm(-vmax, vmax))
    ax.set_xticks(range(len(score_cols)))
    ax.set_xticklabels(score_cols, rotation=90)
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels([f"{i}  (n={int(evidence['n_cells'].iloc[k])})"
                        for k, i in enumerate(matrix.index)])
    ax.set_ylabel("pass-1 cluster")
    # Both flags are marked, because they fail differently. A narrow margin says
    # the panel could not choose between two identities for one population; low
    # support says the cluster is not one population. Showing only the first
    # would repeat the omission that let an 853-cell mixture through unremarked.
    for k, (low, mixed) in enumerate(zip(evidence["low_confidence"],
                                         evidence.get("mixed", False))):
        flag = ("*" if low else "") + ("†" if mixed else "")
        if flag:
            ax.text(len(score_cols) - 0.3, k, flag, va="center", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.12)
    cbar.set_label("marker set score\n(red = enriched, blue = depleted)")
    ax.set_title("Compartment assignment evidence  "
                 "(* narrow margin, † mixed cluster)")
    figs.save_figure(fig_dir, "pass1_compartment_scores")


def plot_composition_points(adata: ad.AnnData, label: str, fig_dir: Path,
                            name: str) -> None:
    """Cell-type proportions per library, one point per library.

    DMM asks for graphs that show the true data spread and for individual data
    points rather than bars. Each point is the proportion of one 10x library;
    with one library per timepoint that is the entire dataset and drawing an
    error bar would imply replication that does not exist, so the count each
    proportion rests on is annotated instead. With several libraries per
    timepoint the points themselves are the spread.
    """
    counts = pd.crosstab(adata.obs["cell_type"], adata.obs["sample"])
    comp = counts / counts.sum(axis=0)
    sample_to_tp = {s: tp for s, tp in _sample_timepoints(adata).items()
                    if s in comp.columns}
    timepoints = list(dict.fromkeys(sample_to_tp.values()))
    per_tp = min(
        sum(1 for t in sample_to_tp.values() if t == tp) for tp in timepoints
    )

    # Long labels are what stops this fitting a single column at 8 pt, so the
    # low-confidence suffix is carried as an asterisk and defined in the legend
    # text instead of being spelled out on every affected row.
    n_low = sum("low confidence" in str(t) for t in comp.index)

    def tick_label(k: int, cell_type: object) -> str:
        """Wrap rather than truncate: cell-type names are the payload of this
        axis, and an elided label is worse than a two-line one."""
        stem = str(cell_type).replace(" (low confidence)", "*")
        wrapped = textwrap.fill(stem.replace("_", " "), width=20)
        return f"{wrapped}\n({int(counts.iloc[k].sum())} cells)"

    labels = [tick_label(k, t) for k, t in enumerate(comp.index)]
    # Row pitch sized to the tallest label: a three-line label needs ~11 mm at
    # 8 pt, and rows closer than their own labels overlap - which is exactly
    # what a fixed per-row height produced.
    max_lines = max(lbl.count("\n") + 1 for lbl in labels)
    pitch = max(0.30, 0.13 * max_lines + 0.07)

    _fig, ax = figs.figure(
        width=figs.SINGLE_COL,
        height=min(0.6 + pitch * len(comp), figs.WORKING_HEIGHT),
    )
    y = np.arange(len(comp))
    markers = ["o", "^", "s", "D"]          # DMM's preferred symbol set
    for j, tp in enumerate(timepoints):
        libs = [s for s, t in sample_to_tp.items() if t == tp]
        band = (j - (len(timepoints) - 1) / 2) * 0.18
        for k, lib in enumerate(libs):
            within = 0.0 if len(libs) == 1 else (k / (len(libs) - 1) - 0.5) * 0.14
            ax.scatter(comp[lib], y + band + within, s=18,
                       marker=markers[j % len(markers)],
                       color=figs.CATEGORICAL[j % len(figs.CATEGORICAL)],
                       label=tp if k == 0 else None,
                       edgecolor="black", linewidth=0.3, zorder=3)
    for k in y:
        ax.hlines(k, 0, float(comp.iloc[k].max()), color="0.85",
                  linewidth=0.5, zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Proportion of library")
    ax.set_xlim(left=0)
    ax.invert_yaxis()
    ax.set_title(f"{label} composition by timepoint", loc="left")
    # Below the x-axis label: the data area can have points at any height, so
    # an in-axes legend has nowhere reliable to go, and next to the title it
    # collides whenever the title fills the single-column width. The tick
    # labels and x-label occupy a fixed ~0.45" (they are set in points), so the
    # anchor offset is converted from inches to axes fraction - a fixed
    # fraction lands far below on tall figures and on top of the label on
    # short ones.
    axes_height = _fig.get_size_inches()[1] * ax.get_position().height
    ax.legend(title=None, loc="upper center",
              bbox_to_anchor=(0.5, -0.45 / axes_height),
              ncols=len(timepoints), frameon=False, borderaxespad=0,
              handletextpad=0.2, columnspacing=1.2)
    # Saved at the width it was designed at. save_figure defaults to a double
    # column, which would silently scale this single-column panel up by 2x and
    # halve its effective type size relative to the page.
    figs.save_figure(fig_dir, name, width=figs.SINGLE_COL)

    replication = (
        f"Each timepoint contributed {per_tp if per_tp > 1 else 'one'} "
        f"librar{'ies' if per_tp > 1 else 'y'}; each point is one library. "
        if per_tp > 1 else
        "One 10x library per timepoint, so there is no biological replication "
        "and no error bars are shown: these proportions are descriptive and "
        "are not tested for a difference between stages. "
    )
    (fig_dir / f"{name}_legend.txt").write_text(
        f"{label} composition by timepoint. Points show the proportion of each "
        f"library assigned to a given cell type; the number in parentheses after "
        f"each label is the total cells of that type across all libraries. "
        + (f"Asterisks mark the {n_low} cell type(s) whose top two marker-set "
           f"scores differed by less than the confidence margin, so the "
           f"assignment could not be resolved by the marker panel alone. "
           if n_low else "")
        + replication + "\n"
    )


def plot_pass2(adata: ad.AnnData, label: str, fig_dir: Path,
               genes: tuple[str, ...] = ()) -> None:
    """Figure set for one re-clustered compartment."""
    _embedding_overview(adata, "umap", fig_dir, f"pass2_{label}_umap")
    _embedding_side_legend(adata, "umap", "cell_type", fig_dir,
                           f"pass2_{label}_cell_types")
    if "X_tsne" in adata.obsm:
        _embedding_overview(adata, "tsne", fig_dir, f"pass2_{label}_tsne")

    _de_dotplot(adata, fig_dir, f"pass2_{label}_dotplot_top_de")

    plot_composition_points(adata, label, fig_dir, f"pass2_{label}_composition")

    if genes:
        # Only plot genes with at least one non-zero count: an all-zero gene is
        # rendered on an artificial symmetric colour scale that looks like real
        # variation. See plot_genes_of_interest for the same guard.
        present = [
            g for g in genes
            if g in adata.raw.var_names and adata.raw[:, g].X.sum() > 0
        ]
        if present:
            figs.scanpy_panel_size(n_cols=3)
            sc.pl.umap(adata, color=present, use_raw=True, ncols=3,
                       cmap=figs.SEQUENTIAL, vmin=0, vmax="p99", frameon=False,
                       wspace=PANEL_SPACING, show=False)
            figs.save_figure(fig_dir, f"pass2_{label}_goi_umap")
            _marker_dotplot(adata, present, "cell_type", fig_dir,
                            f"pass2_{label}_goi_dotplot")


def run(cfg: Config) -> ad.AnnData:
    """Execute the full analysis and write every artefact to `cfg.outdir`.

    Set-up only; the analysis itself lives in `_run` so the whole body executes
    inside `headless_plotting()` without an extra indentation level.
    """
    sc.settings.verbosity = 1
    sc.settings.n_jobs = -1
    # Every sc.pl.* call below passes show=False; this is the belt to that
    # braces, so a plotting helper that forgets cannot re-trigger plt.show().
    sc.settings.autoshow = False
    figs.apply_style()

    cfg.validate_inputs()
    fig_dir = cfg.outdir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = fig_dir

    with headless_plotting():
        return _run(cfg, fig_dir)


def _run(cfg: Config, fig_dir: Path) -> ad.AnnData:
    """Body of `run`, executed inside the headless-plotting context."""
    # Imported here rather than at module scope: pipeline and cluster both
    # reference each other's helpers, and a top-level import would be circular.
    from .cluster import load_markers

    log.info("Stage 1/3: load and QC")
    n_raw: dict[str, int | None] = {s.name: None for s in cfg.samples}
    adata = build_qc_dataset(cfg, raw_counts=n_raw)
    plot_qc(adata, fig_dir)

    log.info("Stage 2/3: normalise, embed, cluster, annotate")
    adata, scores, de = process(adata, cfg)

    log.info("Stage 3/5: pass-1 figures")
    plot_embeddings(adata, fig_dir, cfg.cluster.compute_tsne)
    plot_markers(adata, load_markers(cfg.compartment_markers), fig_dir)

    goi = plot_genes_of_interest(adata, cfg.genes_of_interest, fig_dir)
    if not goi.empty:
        goi.to_csv(cfg.outdir / "genes_of_interest_expression.csv", index=False)

    log.info("Stage 4/5: compartment assignment")
    evidence = assign_compartments(adata, cfg)
    evidence.to_csv(cfg.outdir / "pass1_compartment_evidence.csv")
    plot_compartments(adata, evidence, fig_dir)

    log.info("Stage 5/5: pass-2 sub-clustering")
    subsets: dict[str, ad.AnnData] = {}
    for group in cfg.compartment.subset_groups:
        sub = subset_compartment(adata, group, cfg.compartment.min_cells_for_subset)
        if sub is None:
            continue
        # The ocular compartment gets the fine retinal panel; every other
        # compartment is re-scored against the broad panel, which is the honest
        # ceiling without a curated panel for that tissue.
        panel = cfg.markers_file if group == "ocular" else cfg.compartment_markers
        sub, sub_scores, sub_de = subcluster(sub, cfg, panel, label=group)

        sub_dir = cfg.outdir / group
        sub_fig = sub_dir / "figures"
        sub_dir.mkdir(parents=True, exist_ok=True)
        plot_pass2(sub, group, sub_fig, cfg.genes_of_interest)
        sub_scores.to_csv(sub_dir / "subcluster_celltype_scores.csv")
        sub_de.to_csv(sub_dir / "subcluster_markers_wilcoxon.csv", index=False)
        pd.crosstab(sub.obs["cell_type"], sub.obs["sample"]).to_csv(
            sub_dir / "celltype_counts.csv")
        sub.write_h5ad(sub_dir / f"{group}.h5ad", compression="gzip")
        subsets[group] = sub
        log.info("  [%s] written to %s", group, sub_dir)

        if group == "ocular" and cfg.article_panel:
            # Imported lazily: article builds on this module's helpers, so a
            # top-level import would be circular.
            from .article import plot_article_panels
            plot_article_panels(sub, cfg.article_panel, cfg.markers_file,
                                cfg.outdir / "article")

    scores.to_csv(cfg.outdir / "cluster_celltype_scores.csv")
    de.to_csv(cfg.outdir / "cluster_markers_wilcoxon.csv", index=False)
    adata.write_h5ad(cfg.outdir / "processed.h5ad", compression="gzip")

    summary = (
        adata.obs.groupby(["cell_type", "sample"], observed=True)
        .size().unstack(fill_value=0)
    )
    summary.to_csv(cfg.outdir / "celltype_counts.csv")

    write_manifest(cfg, cfg.outdir / "run_manifest.json", extra={
        "cells_loaded": n_raw,
        "cells_retained": int(adata.n_obs),
        # Recorded explicitly: a run in which Scrublet failed removes exactly as
        # many cells as one in which it found nothing, and the two are otherwise
        # indistinguishable from the outputs.
        "doublet_detection": adata.uns.get("doublet_detection", "unknown"),
        "genes_retained": int(adata.n_vars),
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "cell_types": sorted(map(str, adata.obs["cell_type"].unique())),
        "compartments": adata.obs["compartment_group"].value_counts().to_dict(),
        "pass2_subsets": {
            name: {"cells": int(s.n_obs),
                   "subclusters": int(s.obs["leiden"].nunique()),
                   "cell_types": sorted(map(str, s.obs["cell_type"].unique()))}
            for name, s in subsets.items()
        },
    })

    log.info("Done. Outputs in %s", cfg.outdir.resolve())
    return adata
