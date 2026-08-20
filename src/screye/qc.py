"""Loading and quality control for CellRanger 10x outputs."""

from __future__ import annotations

import logging
import warnings

import anndata as ad
import numpy as np
import scanpy as sc

from .config import Config, QCConfig, SampleSpec

log = logging.getLogger(__name__)


def load_sample(spec: SampleSpec) -> ad.AnnData:
    """Read one CellRanger `filtered_feature_bc_matrix.h5` into an AnnData.

    CellBender's `*_cellbender_filtered.h5` outputs are written in the same
    CellRanger v3 HDF5 layout and load through the same reader, so the
    ambient-corrected per-library matrices published by the upstream Snakemake
    workflow need no separate code path.

    Gene symbols are made unique, and sample-level metadata is attached to
    `.obs` so that provenance survives concatenation.
    """
    log.info("Reading %s from %s", spec.name, spec.h5)

    with warnings.catch_warnings():
        # scanpy warns about duplicate gene symbols while reading; they are made
        # unique on the very next line, so the warning is noise here. The filter
        # is scoped to this one call and matched on the exact message and
        # category, so any other UserWarning from the reader still surfaces.
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = sc.read_10x_h5(spec.h5)

    # Record what was suppressed rather than silently discarding it. Duplicate
    # symbols are normal in Ensembl-based references (several gene IDs mapping
    # to one name), but an unexpectedly large count points at a malformed or
    # mismatched reference GTF.
    n_duplicated = int(adata.var_names.duplicated().sum())
    adata.var_names_make_unique()
    if n_duplicated:
        log.info("  %d duplicate gene symbol(s) made unique", n_duplicated)
    adata.obs["sample"] = spec.name
    adata.obs["timepoint"] = spec.timepoint
    adata.obs_names = [f"{spec.name}_{bc}" for bc in adata.obs_names]
    # Keep the raw counts addressable after normalisation.
    adata.layers["counts"] = adata.X.copy()
    log.info("  %d cells x %d genes", adata.n_obs, adata.n_vars)
    return adata


def annotate_gene_classes(adata: ad.AnnData, cfg: QCConfig) -> None:
    """Flag mitochondrial and ribosomal genes and compute QC metrics.

    Note the lower-case ``mt-`` prefix: zebrafish gene symbols follow ZFIN
    nomenclature, so the human ``MT-`` convention silently matches nothing and
    yields a mitochondrial fraction of exactly zero for every cell. That failure
    mode is silent, which is why it is asserted against below.
    """
    adata.var["mt"] = adata.var_names.str.lower().str.startswith(cfg.mito_prefix)
    adata.var["ribo"] = adata.var_names.str.lower().str.startswith(cfg.ribo_prefixes)

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True
    )

    n_mito = int(adata.var["mt"].sum())
    if n_mito == 0:
        log.warning(
            "No genes matched the mitochondrial prefix %r. Check the species "
            "gene-symbol convention before trusting the QC filter.",
            cfg.mito_prefix,
        )
    else:
        log.info("  %d mitochondrial genes flagged", n_mito)


def mad_outlier(values: np.ndarray, n_mads: float, upper_only: bool = False) -> np.ndarray:
    """Boolean mask of values more than `n_mads` median absolute deviations out.

    Using MADs rather than fixed cut-offs means each sample is filtered against
    its own distribution, which matters when comparing timepoints with different
    sequencing depth.
    """
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros_like(values, dtype=bool)
    deviation = (values - median) / mad
    if upper_only:
        return deviation > n_mads
    return np.abs(deviation) > n_mads


def flag_low_quality(adata: ad.AnnData, cfg: QCConfig) -> None:
    """Mark low-quality cells in `.obs["qc_fail"]` without removing them yet.

    Flag-then-filter keeps the discarded population inspectable: if an entire
    biological cluster is being removed, that is visible in the QC plots rather
    than invisible in the cell count.
    """
    counts = adata.obs["total_counts"].to_numpy()
    genes = adata.obs["n_genes_by_counts"].to_numpy()
    pct_mt = adata.obs["pct_counts_mt"].to_numpy()

    fail = (
        mad_outlier(np.log1p(counts), cfg.n_mads)
        | mad_outlier(np.log1p(genes), cfg.n_mads)
        | mad_outlier(pct_mt, cfg.n_mads, upper_only=True)
        | (pct_mt > cfg.max_pct_mito)
        | (genes < cfg.min_genes_per_cell)
    )
    adata.obs["qc_fail"] = fail
    log.info(
        "  %d/%d cells flagged as low quality (%.1f%%)",
        fail.sum(), adata.n_obs, 100 * fail.mean(),
    )


def _scrublet_batches(adata: ad.AnnData) -> list[tuple[str, float, np.ndarray]]:
    """Per-batch (name, threshold, simulated scores) from Scrublet's own record.

    scanpy nests these under `uns["scrublet"]["batches"]` when `batch_key` is
    used and stores them flat when it is not; both shapes are returned the same
    way here so callers need not care which path ran.
    """
    uns = adata.uns.get("scrublet", {})
    if "batches" in uns:
        return [
            (name, float(rec["threshold"]), np.asarray(rec["doublet_scores_sim"]))
            for name, rec in uns["batches"].items()
            if "threshold" in rec and "doublet_scores_sim" in rec
        ]
    if "threshold" in uns and "doublet_scores_sim" in uns:
        return [("all", float(uns["threshold"]),
                 np.asarray(uns["doublet_scores_sim"]))]
    return []


def _apply_threshold(scored: ad.AnnData, cfg: QCConfig) -> str:
    """Validate Scrublet's automatic threshold, or replace it with an explicit one.

    Returns a short note for the manifest. Reports *sensitivity* — the fraction
    of Scrublet's own simulated doublets that the threshold catches — because a
    call rate on its own reads as a doublet rate, and here it is not one: at a
    threshold of 0.25 only about a sixth of simulated doublets are detectable,
    so the cells removed are the detectable minority rather than all of them.
    """
    batches = _scrublet_batches(scored)
    setting = cfg.doublet_threshold

    if isinstance(setting, str) and setting.strip().lower() == "auto":
        # A threshold is degenerate when it sits above every observed cell: no
        # barcode can be called, so "0 doublets" describes the threshold rather
        # than the data. Testing the observed scores is what matters — testing
        # only the simulated ones is not enough, because Scrublet's automatic
        # threshold can leave a handful of extreme simulants above it (12 of
        # 102,462 on GSE158142) while still being unreachable by any real cell.
        observed = scored.obs["doublet_score"].to_numpy()
        highest = float(np.nanmax(observed)) if observed.size else 0.0
        degenerate = [
            f"{name}: threshold {thr:.3f} exceeds every observed cell "
            f"(highest {highest:.3f}); it catches only "
            f"{100 * (sim > thr).mean():.2f}% of simulated doublets"
            for name, thr, sim in batches
            if sim.size and thr >= highest
        ]
        if degenerate:
            raise ValueError(
                "Scrublet's automatic threshold is degenerate — it excludes every "
                "simulated doublet, so no cell can be called:\n  "
                + "\n  ".join(degenerate)
                + "\n\nThis happens when the simulated-doublet score histogram is "
                "unimodal (shallow data, or a tissue that is a continuum of "
                "related states) and Otsu's method has no trough to find. Set "
                "qc.doublet_threshold to an explicit value chosen from that "
                "histogram, or set qc.run_doublet_detection to false and say so "
                "in the write-up."
            )
        thresholds = ", ".join(f"{name}={thr:.3f}" for name, thr, _ in batches)
        return f"threshold auto ({thresholds})"

    threshold = float(setting)
    scored.obs["predicted_doublet"] = (
        scored.obs["doublet_score"].to_numpy() > threshold
    )
    if batches:
        pooled = np.concatenate([sim for _, _, sim in batches if sim.size])
        sensitivity = float((pooled > threshold).mean())
        log.info("  threshold %.3f fixed by config; detects %.1f%% of simulated doublets",
                 threshold, 100 * sensitivity)
        return f"threshold {threshold} (explicit); sensitivity {sensitivity:.3f}"
    return f"threshold {threshold} (explicit)"


def detect_doublets(adata: ad.AnnData, cfg: QCConfig, seed: int) -> str:
    """Run Scrublet per sample and record scores and calls.

    Doublets are scored within each sample: pooling samples first would let
    barcodes from one library act as simulated doublet partners for another.

    Scoring is restricted to cells that pass QC. That is partly a correctness
    choice and partly a hard requirement:

    - Simulated doublets are built by summing pairs of observed profiles.
      Including barcodes that are about to be discarded means synthesising
      doublets out of noise and calibrating the threshold against them.
    - Ambient-corrected input makes it mandatory. CellBender can strip a
      barcode down to almost nothing — in the 22-library GSE158142 set, 41
      barcodes end at exactly zero counts and 178 at two or fewer. Scrublet's
      own preprocessing drops anything below three counts, after which scanpy
      cannot reindex the per-batch results back onto the full set and raises
      `KeyError: [...] not in index`. Cell Ranger's cell calls never contain
      such barcodes, so this failure appears only after ambient correction.

    `min_genes_per_cell` (200 by default) puts every cell reaching Scrublet far
    above its internal threshold.

    Returns a short status string, recorded in the run manifest so that a run
    which skipped doublet detection cannot be mistaken for one that found none.
    """
    adata.obs["predicted_doublet"] = False
    adata.obs["doublet_score"] = np.nan

    if not cfg.run_doublet_detection:
        return "disabled"

    passing = ~adata.obs["qc_fail"].to_numpy()
    if passing.sum() == 0:
        log.warning("  no cells pass QC; skipping doublet detection")
        return "skipped: no cells passed QC"

    scored = adata[passing].copy()
    log.info("  scoring %d QC-passing cells (%d excluded)",
             scored.n_obs, adata.n_obs - scored.n_obs)

    try:
        sc.pp.scrublet(
            scored,
            batch_key="sample" if scored.obs["sample"].nunique() > 1 else None,
            expected_doublet_rate=cfg.expected_doublet_rate,
            random_state=seed,
        )
    except Exception as exc:  # noqa: BLE001 # pragma: no cover - scrublet fails on tiny inputs
        # Deliberately not fatal: Scrublet legitimately cannot converge on very
        # small inputs, and the synthetic test fixtures are one such case. But
        # the outcome is returned rather than swallowed, because a silent
        # no-op here produces a complete-looking result set in which no doublet
        # was ever removed.
        log.error("Scrublet failed (%s); CONTINUING WITHOUT DOUBLET CALLS", exc)
        return f"failed: {type(exc).__name__}: {str(exc)[:200]}"

    # Outside the guard above, deliberately. A degenerate or mis-set threshold
    # is a configuration error the caller must resolve, not a convergence
    # failure to be tolerated — swallowing it would reproduce exactly the silent
    # zero-doublet run this function exists to prevent.
    note = _apply_threshold(scored, cfg)

    adata.obs.loc[scored.obs_names, "doublet_score"] = scored.obs["doublet_score"]
    adata.obs.loc[scored.obs_names, "predicted_doublet"] = (
        scored.obs["predicted_doublet"].astype(bool)
    )
    n_doublets = int(adata.obs["predicted_doublet"].sum())
    log.info("  %d cells called as doublets (%.2f%% of QC-passing)",
             n_doublets, 100 * n_doublets / scored.n_obs)
    return f"ok: {n_doublets} called on {scored.n_obs} QC-passing cells; {note}"


def apply_filters(adata: ad.AnnData, cfg: QCConfig) -> ad.AnnData:
    """Remove flagged cells and rarely-detected genes."""
    before = adata.n_obs
    keep = ~adata.obs["qc_fail"].to_numpy() & ~adata.obs["predicted_doublet"].to_numpy()
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=cfg.min_cells_per_gene)
    log.info("  retained %d/%d cells, %d genes", adata.n_obs, before, adata.n_vars)
    return adata


def build_qc_dataset(cfg: Config, raw_counts: dict[str, int] | None = None) -> ad.AnnData:
    """Load every sample, QC it, and return the concatenated object.

    `raw_counts`, if given, is populated with the pre-filter cell count per
    sample so that attrition can be reported in the run manifest.
    """
    adatas = [load_sample(spec) for spec in cfg.samples]
    if raw_counts is not None:
        for spec, a in zip(cfg.samples, adatas):
            raw_counts[spec.name] = int(a.n_obs)
    adata = ad.concat(adatas, label=None, index_unique=None, join="inner")
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    adata.obs["timepoint"] = adata.obs["timepoint"].astype("category")
    log.info("Concatenated: %d cells x %d genes", adata.n_obs, adata.n_vars)

    annotate_gene_classes(adata, cfg.qc)
    flag_low_quality(adata, cfg.qc)
    # Kept in `.uns` so the status travels with the object into processed.h5ad
    # as well as into the run manifest.
    adata.uns["doublet_detection"] = detect_doublets(adata, cfg.qc, cfg.seed)
    return apply_filters(adata, cfg.qc)
