"""Compartment assignment and subsetting for two-pass annotation.

Pass 1 clusters the whole dissociation and assigns each cluster to a broad tissue
compartment. Pass 2 subsets one compartment and re-runs feature selection,
embedding and clustering within it.

Why two passes
--------------
Highly variable gene selection is a global operation. On a whole-head
dissociation the top 2,000 HVGs are dominated by the axes that separate brain
from muscle from blood, and the genes that distinguish a bipolar cell from an
amacrine cell are not among them. Subsetting and re-selecting recovers that
structure. The same argument appears in the focused-clustering analyses of Sur
et al. 2023 (Dev Cell 58:3028) and in Raj et al. 2020 (Neuron 108:1058), the
source of these data.

The cost is that compartment assignment is now a decision the rest of the
analysis depends on, so it is made explicitly, recorded, and reported with the
evidence behind it rather than being buried inside a labelling step.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from .config import Config

log = logging.getLogger(__name__)


def load_panel(path: str | Path) -> dict[str, list[str]]:
    """Read a marker panel from YAML."""
    with Path(path).open() as handle:
        return yaml.safe_load(handle)["cell_types"]


def score_panel(
    adata: ad.AnnData,
    panel: dict[str, list[str]],
    groupby: str,
    prefix: str = "score",
) -> pd.DataFrame:
    """Score every cell against each marker set; return the cluster mean matrix.

    `sc.tl.score_genes` compares mean expression of the set against a randomly
    sampled reference set matched for expression bin, so scores are comparable
    between sets of different size and expression level. They are NOT
    probabilities and should not be read as such.
    """
    present = {
        name: [g for g in genes if g in adata.raw.var_names]
        for name, genes in panel.items()
    }
    empty = [name for name, genes in present.items() if not genes]
    if empty:
        log.warning(
            "No markers found in the data for: %s. Check symbols against ZFIN "
            "and against the reference GTF.", ", ".join(empty),
        )

    for name, genes in present.items():
        if genes:
            sc.tl.score_genes(adata, genes, score_name=f"{prefix}_{name}", use_raw=True)

    cols = [c for c in adata.obs.columns if c.startswith(f"{prefix}_")]
    return (
        adata.obs.groupby(groupby, observed=True)[cols].mean()
        .rename(columns=lambda c: c.removeprefix(f"{prefix}_"))
    )


def assign_by_margin(
    scores: pd.DataFrame,
    min_margin: float = 0.05,
) -> tuple[pd.Series, pd.Series]:
    """Label each row by its top-scoring column, flagging narrow decisions.

    Returns (labels, margins). A margin is the difference between the best and
    second-best score. A small margin means the panel cannot separate those two
    identities for that cluster; asserting the winner would be overconfident, so
    the label is suffixed rather than silently accepted.
    """
    best = scores.idxmax(axis=1)
    margins = scores.apply(
        lambda row: float(np.diff(np.sort(row.to_numpy())[-2:])[0])
        if len(row) > 1 else np.inf,
        axis=1,
    )
    labels = pd.Series(
        {
            idx: (lab if margins[idx] >= min_margin else f"{lab} (low confidence)")
            for idx, lab in best.items()
        },
        name="label",
    )
    return labels, margins


def assign_compartments(
    adata: ad.AnnData,
    cfg: Config,
    groupby: str = "leiden",
) -> pd.DataFrame:
    """Assign each pass-1 cluster to a tissue compartment.

    Writes `.obs["compartment"]` (fine label, e.g. "ocular_photoreceptor") and
    `.obs["compartment_group"]` (coarse, e.g. "ocular"). Returns the evidence
    table: score per compartment per cluster, plus the margin and cell count.
    """
    panel = load_panel(cfg.compartment_markers)
    scores = score_panel(adata, panel, groupby=groupby, prefix="cscore")
    labels, margins = assign_by_margin(scores, cfg.compartment.min_margin)

    adata.obs["compartment"] = adata.obs[groupby].map(labels).astype("category")

    # Coarse grouping. Fine labels carry the compartment as their prefix, so the
    # mapping is mechanical rather than another hand-curated list to drift.
    def coarse(label: str) -> str:
        stem = label.replace(" (low confidence)", "")
        if stem.startswith("ocular_"):
            return "ocular"
        if stem in {"neuron_differentiated", "neural_progenitor_radial_glia",
                    "oligodendrocyte_opc"}:
            return "neural"
        return "other"

    adata.obs["compartment_group"] = (
        adata.obs["compartment"].astype(str).map(coarse).astype("category")
    )

    evidence = scores.copy()
    evidence["assigned"] = labels
    evidence["margin"] = margins
    evidence["n_cells"] = adata.obs[groupby].value_counts()
    evidence["low_confidence"] = evidence["assigned"].str.contains("low confidence")

    n_low = int(evidence["low_confidence"].sum())
    counts = adata.obs["compartment_group"].value_counts().to_dict()
    log.info("  compartment assignment: %s", counts)
    if n_low:
        log.warning(
            "  %d/%d cluster(s) assigned on a margin below %.2f; treat those "
            "compartment calls as provisional",
            n_low, len(evidence), cfg.compartment.min_margin,
        )
    return evidence.sort_values("margin")


def subset_compartment(
    adata: ad.AnnData,
    group: str,
    min_cells: int = 200,
) -> ad.AnnData | None:
    """Return the cells of one compartment, restored to pre-HVG state.

    The returned object carries the full log-normalised gene set from `.raw`, not
    the whole-tissue HVG subset. Re-running feature selection on the subset is
    the entire point of the second pass, so it must start from all genes.
    """
    mask = adata.obs["compartment_group"] == group
    n = int(mask.sum())
    if n < min_cells:
        log.warning(
            "Compartment %r has %d cells (< %d); skipping second pass. "
            "Sub-clustering below this size produces clusters that reflect "
            "sampling noise rather than biology.", group, n, min_cells,
        )
        return None

    subset = adata.raw.to_adata()[mask].copy()
    # Carry forward the pass-1 assignment so pass-2 clusters can be compared
    # against it, and drop the whole-tissue scores which no longer apply.
    for key in ("sample", "timepoint", "leiden", "compartment", "compartment_group"):
        if key in adata.obs:
            subset.obs[key + ("_pass1" if key in ("leiden",) else "")] = adata.obs[key][mask]
    subset.obs = subset.obs.loc[:, ~subset.obs.columns.str.startswith(("cscore_", "score_"))]

    # Drop pass-1 derived state. A dendrogram, colour map or rank_genes_groups
    # result computed on the whole-tissue clustering is not merely stale here,
    # it is wrong: scanpy will reuse it against the new cluster labels and
    # either error or, worse, silently plot pass-1 groupings under pass-2 names.
    for key in ("dendrogram_leiden", "rank_genes_groups", "leiden",
                "leiden_colors", "cell_type_colors", "compartment_colors",
                "compartment_group_colors", "neighbors", "umap", "pca", "tsne",
                "log1p", "hvg"):
        subset.uns.pop(key, None)
    subset.obsm.clear()
    subset.varm.clear()
    subset.obsp.clear()

    log.info("  %s compartment: %d cells x %d genes carried into pass 2",
             group, subset.n_obs, subset.n_vars)
    return subset
