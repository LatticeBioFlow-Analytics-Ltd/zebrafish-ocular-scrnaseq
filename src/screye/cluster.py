"""Normalisation, dimensionality reduction, clustering and annotation."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import yaml

from .config import ClusterConfig, Config

log = logging.getLogger(__name__)


def normalise(adata: ad.AnnData, cfg: ClusterConfig) -> ad.AnnData:
    """Library-size normalise, log-transform, and select highly variable genes.

    The full log-normalised matrix is kept in `.raw` so that marker expression
    can be plotted for genes excluded from the HVG set.
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    batch_key = cfg.batch_key if adata.obs[cfg.batch_key].nunique() > 1 else None
    sc.pp.highly_variable_genes(
        adata, n_top_genes=cfg.n_top_genes, batch_key=batch_key, flavor="seurat"
    )
    log.info("  %d highly variable genes", int(adata.var["highly_variable"].sum()))
    return adata


def embed(adata: ad.AnnData, cfg: ClusterConfig, seed: int) -> ad.AnnData:
    """PCA, optional batch correction, neighbour graph, UMAP and t-SNE.

    Batch correction is deliberately opt-in. When the samples are different
    developmental stages, `sample` and `stage` are the same variable, so
    integrating on `sample` removes the developmental signal being measured.
    Correct only if you have a genuine technical batch to remove.
    """
    # Zero-centring cannot preserve sparsity, so this materialises a dense
    # array. Reported rather than suppressed silently: at large cell counts this
    # is the step that will exhaust memory, and knowing the size in advance is
    # more useful than discovering it as a kill signal.
    log.info("  scaling %d x %d (~%.1f GB dense)", adata.n_obs, adata.n_vars,
             adata.n_obs * adata.n_vars * 4 / 1e9)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="zero-centering a sparse array",
                                category=UserWarning)
        sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=cfg.n_pcs, svd_solver="arpack", random_state=seed)

    use_rep = "X_pca"
    if cfg.batch_correct and adata.obs[cfg.batch_key].nunique() > 1:
        log.warning(
            "Batch-correcting on %r. Confirm this is a technical batch and not "
            "the biological contrast of interest.", cfg.batch_key,
        )
        sc.external.pp.harmony_integrate(adata, key=cfg.batch_key)
        use_rep = "X_pca_harmony"

    sc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, n_pcs=cfg.n_pcs,
                    use_rep=use_rep, random_state=seed)
    sc.tl.umap(adata, random_state=seed)

    if cfg.compute_tsne:
        # UMAP and t-SNE are both reported because they disagree usefully:
        # t-SNE tends to fragment continuous trajectories that UMAP joins.
        # Neither preserves global inter-cluster distance, so cluster
        # adjacency in either embedding is not evidence of relatedness.
        sc.tl.tsne(adata, n_pcs=cfg.n_pcs, use_rep=use_rep,
                   perplexity=cfg.tsne_perplexity, random_state=seed)
    return adata


def cluster(adata: ad.AnnData, cfg: ClusterConfig, seed: int) -> ad.AnnData:
    """Leiden clustering across a resolution sweep.

    Several resolutions are stored so the choice can be justified from the data
    rather than asserted; `leiden` holds the primary resolution.
    """
    for res in cfg.leiden_resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key,
                     random_state=seed, flavor="igraph", n_iterations=2,
                     directed=False)
        log.info("  resolution %.2f -> %d clusters", res, adata.obs[key].nunique())

    adata.obs["leiden"] = adata.obs[f"leiden_{cfg.primary_resolution}"]
    return adata


def load_markers(path: Path) -> dict[str, list[str]]:
    """Read the curated marker panel."""
    with Path(path).open() as handle:
        return yaml.safe_load(handle)["cell_types"]


def score_cell_types(adata: ad.AnnData, markers: dict[str, list[str]]) -> pd.DataFrame:
    """Score each cell against every marker set and label clusters by mean score.

    Returns the cluster-by-cell-type score matrix so that the assignment is
    auditable: a cluster whose top two scores are close is a cluster whose label
    should not be trusted without further evidence.
    """
    present = {
        name: [g for g in genes if g in adata.raw.var_names]
        for name, genes in markers.items()
    }
    for name, genes in present.items():
        if not genes:
            log.warning("No markers for %r found in the data; skipping", name)
            continue
        sc.tl.score_genes(adata, genes, score_name=f"score_{name}", use_raw=True)

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    scores = (
        adata.obs.groupby("leiden", observed=True)[score_cols].mean()
        .rename(columns=lambda c: c.removeprefix("score_"))
    )

    best = scores.idxmax(axis=1)
    top2 = scores.apply(lambda r: r.nlargest(2).to_numpy(), axis=1)
    margin = top2.apply(lambda v: v[0] - v[1] if len(v) > 1 else float("inf"))

    label_map = {
        cl: (lab if margin[cl] > 0.05 else f"{lab} (low confidence)")
        for cl, lab in best.items()
    }
    adata.obs["cell_type"] = adata.obs["leiden"].map(label_map).astype("category")

    n_low = sum("low confidence" in v for v in label_map.values())
    if n_low:
        log.warning("%d cluster(s) labelled with a narrow score margin", n_low)
    return scores


def rank_markers(adata: ad.AnnData) -> pd.DataFrame:
    """Wilcoxon differential expression per cluster, returned as a tidy frame."""
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    return sc.get.rank_genes_groups_df(adata, group=None)


def process(adata: ad.AnnData, cfg: Config) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """Run the full downstream analysis."""
    adata = normalise(adata, cfg.cluster)
    hvg = adata[:, adata.var["highly_variable"]].copy()
    # `.raw` yields a Raw object; converting back to AnnData keeps the full
    # log-normalised matrix addressable for marker plotting after subsetting.
    hvg.raw = adata.raw.to_adata()
    hvg = embed(hvg, cfg.cluster, cfg.seed)
    hvg = cluster(hvg, cfg.cluster, cfg.seed)
    scores = score_cell_types(hvg, load_markers(cfg.markers_file))
    de = rank_markers(hvg)
    return hvg, scores, de


def subcluster(
    adata: ad.AnnData,
    cfg: Config,
    markers_file: Path,
    label: str,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """Pass 2: re-run feature selection, embedding and clustering on a subset.

    Deliberately re-does HVG selection rather than reusing the whole-tissue set.
    Feature selection is global, so on a whole-head dissociation the top genes
    separate brain from muscle from blood, and the genes distinguishing a
    bipolar cell from an amacrine cell do not make the list. Re-selecting within
    the compartment is what recovers that structure, and is the entire reason
    for the second pass.

    `label` names the compartment and is used only for logging.
    """
    from .compartments import assign_by_margin, load_panel, score_panel

    ccfg = cfg.compartment
    log.info("  [%s] re-selecting features on %d cells", label, adata.n_obs)

    # NOT re-normalised. The subset arrives already log-normalised from pass 1,
    # and library-size normalisation is a per-cell operation, so subsetting
    # cells does not invalidate it. Re-running normalize_total/log1p here would
    # log-transform twice and silently compress the expression range.
    adata.raw = adata

    batch_key = "sample" if adata.obs["sample"].nunique() > 1 else None
    sc.pp.highly_variable_genes(
        adata, n_top_genes=min(ccfg.n_top_genes, adata.n_vars - 1),
        batch_key=batch_key, flavor="seurat",
    )
    hvg = adata[:, adata.var["highly_variable"]].copy()
    hvg.raw = adata.raw.to_adata()
    log.info("  [%s] %d highly variable genes within compartment",
             label, int(adata.var["highly_variable"].sum()))

    n_pcs = min(ccfg.n_pcs, hvg.n_obs - 1, hvg.n_vars - 1)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="zero-centering a sparse array",
                                category=UserWarning)
        sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=n_pcs, svd_solver="arpack", random_state=cfg.seed)
    sc.pp.neighbors(hvg, n_neighbors=min(ccfg.n_neighbors, hvg.n_obs - 1),
                    n_pcs=n_pcs, random_state=cfg.seed)
    sc.tl.umap(hvg, random_state=cfg.seed)
    if cfg.cluster.compute_tsne:
        sc.tl.tsne(hvg, n_pcs=n_pcs,
                   perplexity=min(cfg.cluster.tsne_perplexity, (hvg.n_obs - 1) / 3),
                   random_state=cfg.seed)

    for res in ccfg.leiden_resolutions:
        sc.tl.leiden(hvg, resolution=res, key_added=f"leiden_{res}",
                     random_state=cfg.seed, flavor="igraph", n_iterations=2,
                     directed=False)
    hvg.obs["leiden"] = hvg.obs[f"leiden_{ccfg.primary_resolution}"]
    log.info("  [%s] %d sub-clusters at resolution %.2f",
             label, hvg.obs["leiden"].nunique(), ccfg.primary_resolution)

    panel = load_panel(markers_file)
    scores = score_panel(hvg, panel, groupby="leiden", prefix="score")
    labels, margins = assign_by_margin(scores, ccfg.min_margin)
    hvg.obs["cell_type"] = hvg.obs["leiden"].map(labels).astype("category")

    scores = scores.copy()
    scores["assigned"] = labels
    scores["margin"] = margins
    scores["n_cells"] = hvg.obs["leiden"].value_counts()

    n_low = int(labels.str.contains("low confidence").sum())
    if n_low:
        log.warning("  [%s] %d/%d sub-cluster(s) labelled on a narrow margin",
                    label, n_low, len(labels))

    de = rank_markers(hvg)
    return hvg, scores, de
