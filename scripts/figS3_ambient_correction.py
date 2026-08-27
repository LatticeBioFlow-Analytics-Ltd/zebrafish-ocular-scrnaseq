"""Fig. S3 panel: what ambient-RNA correction did to the RPE cluster.

The other S3 panels describe one run. This one is a comparison, so it cannot be
produced inside the pipeline: it needs the corrected and the uncorrected analysis
side by side, and each takes a full run to make.

The argument it supports is narrow and specific. Detection of the three
transporters in the RPE could in principle be ambient contamination — soup from
abundant transcripts elsewhere in the eye landing in RPE droplets. If it were,
correction would remove it. So the panel plots the same cells before and after,
and asks which transcripts survive:

  A  Detection in matched RPE cells, uncorrected against corrected, per panel
     gene. Ambient sentinels should fall; genuinely expressed transcripts
     should not.
  B  Counts removed per library, so the magnitude of the correction is visible
     and its variation between libraries is not hidden by an average.

Cells are matched by barcode across the two runs. The uncorrected analysis used
per-timepoint aggregates, whose barcodes carry a GEM-group suffix (-1..-N)
assigned by `cellranger aggr` in the row order of its input CSV; that order is
the sorted library list per timepoint. The mapping is verified here rather than
assumed, because getting it wrong would silently compare unrelated cells.

    pixi run -e dev python scripts/figS3_ambient_correction.py
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from screye import figures as figs  # noqa: E402
from screye.article import load_panel_spec  # noqa: E402

warnings.filterwarnings("ignore")
log = logging.getLogger("figS3_ambient")

MIN_OVERLAP = 0.5  # a correct GEM->library mapping overlaps ~0.8-0.99


def libraries_by_stage(sheet: Path) -> dict[str, list[str]]:
    """Sorted library list per timepoint — the order `aggr` assigns GEM groups in."""
    out = defaultdict(list)
    with sheet.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("library"):
                out[row["timepoint"]].append(row["library"])
    return {k: sorted(v) for k, v in out.items()}


def key_corrected(name: str) -> tuple[str, str]:
    library, barcode = name.rsplit("_", 1)
    return library, barcode.split("-")[0]


def key_uncorrected(name: str, libs: dict[str, list[str]],
                    stage_of: dict[str, str]) -> tuple[str, str]:
    sample, rest = name.split("_", 1)
    barcode, gem = rest.rsplit("-", 1)
    return libs[stage_of[sample]][int(gem) - 1], barcode


def verify_mapping(old: ad.AnnData, new: ad.AnnData, libs, stage_of) -> None:
    """Confirm each GEM group matches its assumed library by barcode overlap."""
    by_library = defaultdict(set)
    for n in new.obs_names:
        library, barcode = key_corrected(n)
        by_library[library].add(barcode)

    by_gem = defaultdict(set)
    for n in old.obs_names:
        sample, rest = n.split("_", 1)
        barcode, gem = rest.rsplit("-", 1)
        by_gem[(sample, int(gem))].add(barcode)

    for (sample, gem), barcodes in sorted(by_gem.items()):
        candidates = libs[stage_of[sample]]
        scored = sorted(
            ((len(barcodes & by_library[lib]) / len(barcodes), lib)
             for lib in candidates), reverse=True)
        best_fraction, best_library = scored[0]
        assumed = candidates[gem - 1]
        if best_library != assumed or best_fraction < MIN_OVERLAP:
            raise SystemExit(
                f"GEM-group mapping failed for {sample}-{gem}: best match is "
                f"{best_library} at {best_fraction:.1%}, but position {gem} in "
                f"the sorted library list is {assumed}. Comparing these runs "
                "would pair unrelated cells."
            )
    log.info("GEM-group mapping verified across %d groups", len(by_gem))


def detection(adata: ad.AnnData, genes: list[str], mask: np.ndarray) -> pd.Series:
    """Percent of masked cells with any detected transcript, per gene.

    Read off `.raw`, which holds the log-normalised matrix. Normalisation and
    log1p both map zero to zero and non-zero to non-zero, so detection is
    identical to what it would be on counts — and identical to the definition
    the other S3 panels use.
    """
    out = {}
    for gene in genes:
        if gene not in adata.raw.var_names:
            out[gene] = np.nan
            continue
        column = adata.raw[:, gene].X
        dense = (column.toarray() if hasattr(column, "toarray")
                 else np.asarray(column)).ravel()
        out[gene] = 100 * (dense[mask] > 0).mean()
    return pd.Series(out)


def counts_removed(matched: set, libs, stage_of, data_dir: Path) -> pd.DataFrame:
    """Percent of UMIs removed per library, over cells matched in both runs."""
    before, after = {}, {}
    for sample in ("day5", "day8"):
        path = data_dir / f"{sample}_filtered_feature_bc_matrix.h5"
        if not path.exists():
            log.warning("%s missing; panel B will be omitted", path.name)
            return pd.DataFrame()
        a = sc.read_10x_h5(path)
        a.var_names_make_unique()
        totals = np.asarray(a.X.sum(axis=1)).ravel()
        for name, total in zip(a.obs_names, totals):
            barcode, gem = name.rsplit("-", 1)
            before[(libs[stage_of[sample]][int(gem) - 1], barcode)] = total

    for library in sorted({lib for lib, _ in matched}):
        a = sc.read_10x_h5(data_dir / "cellbender" /
                           f"{library}_cellbender_filtered.h5")
        a.var_names_make_unique()
        totals = np.asarray(a.X.sum(axis=1)).ravel()
        for name, total in zip(a.obs_names, totals):
            after[(library, name.split("-")[0])] = total

    rows = []
    for library, barcode in matched:
        b, a_ = before.get((library, barcode)), after.get((library, barcode))
        if b and a_ is not None and b > 0:
            rows.append({"library": library, "stage": stage_of_library(library),
                         "removed": 100 * (b - a_) / b})
    return pd.DataFrame(rows)


_STAGE_OF_LIBRARY: dict[str, str] = {}


def stage_of_library(library: str) -> str:
    return _STAGE_OF_LIBRARY.get(library, "?")


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)-7s | %(message)s")
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corrected", default="results_cellbender")
    parser.add_argument("--uncorrected", default="results")
    parser.add_argument("--outdir", default="results_cellbender/article/figures")
    args = parser.parse_args()

    libs = libraries_by_stage(root / "config/libraries.tsv")
    for stage, names in libs.items():
        for name in names:
            _STAGE_OF_LIBRARY[name] = stage
    stage_of = {"day5": "5dpf", "day8": "8dpf"}

    new = ad.read_h5ad(root / args.corrected / "ocular/ocular.h5ad")
    old = ad.read_h5ad(root / args.uncorrected / "ocular/ocular.h5ad")
    verify_mapping(old, new, libs, stage_of)

    spec = load_panel_spec(root / "config/article_panel.yaml")
    roles = {g: role for role, genes in spec["panel"].items() for g in genes}
    genes = [g for role in spec["panel"] for g in spec["panel"][role]]

    def rpe_keys(adata, keyfn):
        is_rpe = adata.obs["cell_type"].astype(str).str.startswith(
            spec["rpe_cell_type"]).to_numpy()
        return {keyfn(n) for n, r in zip(adata.obs_names, is_rpe) if r}, is_rpe

    new_rpe, new_mask = rpe_keys(new, key_corrected)
    old_rpe, old_mask = rpe_keys(
        old, lambda n: key_uncorrected(n, libs, stage_of))
    matched = new_rpe & old_rpe
    log.info("RPE cells: %d corrected, %d uncorrected, %d matched",
             len(new_rpe), len(old_rpe), len(matched))
    if len(matched) < 50:
        raise SystemExit("too few matched RPE cells to compare")

    new_sel = np.array([key_corrected(n) in matched for n in new.obs_names])
    old_sel = np.array([key_uncorrected(n, libs, stage_of) in matched
                        for n in old.obs_names])
    table = pd.DataFrame({
        "before": detection(old, genes, old_sel),
        "after": detection(new, genes, new_sel),
    })
    table["role"] = [roles[g] for g in table.index]

    removed = counts_removed(matched, libs, stage_of, root / "data")

    # --- draw -----------------------------------------------------------------
    figs.apply_style()
    width = figs.DOUBLE_COL
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(width, width * 0.42),
        gridspec_kw={"width_ratios": (1.8, 1)})

    order = [g for g in genes if not table.loc[g, ["before", "after"]].isna().all()]
    palette = dict(zip(spec["panel"], figs.CATEGORICAL))
    y = np.arange(len(order))[::-1]
    for yi, gene in zip(y, order):
        b, a_ = table.loc[gene, "before"], table.loc[gene, "after"]
        colour = palette[table.loc[gene, "role"]]
        ax_a.plot([b, a_], [yi, yi], color=colour, lw=1.0, alpha=0.6, zorder=1)
        ax_a.scatter([b], [yi], s=16, facecolor="white", edgecolor=colour,
                     lw=0.9, zorder=2)
        ax_a.scatter([a_], [yi], s=16, color=colour, zorder=3)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(order, style="italic")
    ax_a.set_xlabel("RPE cells detecting the gene (%)")
    ax_a.set_xlim(left=0)
    handles = [plt.Line2D([], [], marker="o", ls="", markerfacecolor=c,
                          markeredgecolor=c, markersize=4, label=r)
               for r, c in palette.items()]
    handles += [
        plt.Line2D([], [], marker="o", ls="", markerfacecolor="white",
                   markeredgecolor="0.4", markersize=4, label="uncorrected"),
        plt.Line2D([], [], marker="o", ls="", color="0.4", markersize=4,
                   label="corrected"),
    ]
    ax_a.legend(handles=handles, frameon=False, fontsize=6,
                loc="lower right", handletextpad=0.2, labelspacing=0.25)
    figs.panel_label(ax_a, "A")

    if not removed.empty:
        summary = (removed.groupby(["stage", "library"])["removed"]
                   .median().reset_index().sort_values(["stage", "removed"]))
        shades = {"5dpf": figs.CATEGORICAL[0], "8dpf": figs.CATEGORICAL[1]}
        for i, (_, r) in enumerate(summary.iterrows()):
            ax_b.scatter(r["removed"], i, s=14, color=shades.get(r["stage"], "0.4"))
        ax_b.set_yticks(range(len(summary)))
        ax_b.set_yticklabels(summary["library"], fontsize=5)
        ax_b.set_xlabel("UMIs removed per cell (%, median)")
        ax_b.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=c,
                                        markersize=4, label=s)
                             for s, c in shades.items()],
                    frameon=False, fontsize=6, loc="lower right")
    else:
        ax_b.set_axis_off()
        ax_b.text(0.5, 0.5, "uncorrected per-library\nmatrices unavailable",
                  ha="center", va="center", fontsize=7, color="0.4")
    figs.panel_label(ax_b, "B")

    out = root / args.outdir
    out.mkdir(parents=True, exist_ok=True)
    figs.save_figure(out, "figS3_ambient_correction", width=width)
    table.round(2).to_csv(out.parent / "figS3_ambient_correction.csv")

    sentinels = spec["panel"]["ambient sentinel"]
    targets = spec["panel"]["target"]
    log.info("ambient sentinels in RPE: %s",
             ", ".join(f"{g} {table.loc[g,'before']:.1f}->{table.loc[g,'after']:.1f}%"
                       for g in sentinels if g in table.index))
    log.info("targets in RPE:           %s",
             ", ".join(f"{g} {table.loc[g,'before']:.1f}->{table.loc[g,'after']:.1f}%"
                       for g in targets if g in table.index))
    return 0


if __name__ == "__main__":
    sys.exit(main())
