"""Unit tests for the QC logic that is easy to get silently wrong."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from screye.config import Config, QCConfig, SampleSpec
from screye.qc import (
    _apply_threshold,
    annotate_gene_classes,
    detect_doublets,
    flag_low_quality,
    mad_outlier,
)


def _toy_adata(n_cells: int = 200, seed: int = 0) -> ad.AnnData:
    """A small but realistically-shaped matrix: mitochondrial genes are a few
    percent of the panel, as in a real library, so percentage-based thresholds
    behave the way they would on real data."""
    rng = np.random.default_rng(seed)
    genes = ["mt-nd1", "mt-co1", "rps3", "rpl7", "rho", "gnat1", "vsx1", "pcna"]
    genes += [f"gene{i:03d}" for i in range(92)]
    counts = rng.negative_binomial(5, 0.3, size=(n_cells, len(genes))).astype(float)
    adata = ad.AnnData(counts)
    adata.var_names = genes
    adata.obs["sample"] = "test"
    return adata


def test_mad_outlier_flags_extremes_symmetrically():
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.normal(10.0, 1.0, 100), [1000.0], [-1000.0]])
    mask = mad_outlier(values, n_mads=5.0)
    assert mask[-1] and mask[-2]
    assert not mask[:100].any()


def test_mad_outlier_upper_only_ignores_low_tail():
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.normal(10.0, 1.0, 100), [1000.0], [-1000.0]])
    mask = mad_outlier(values, n_mads=5.0, upper_only=True)
    assert mask[-2] and not mask[-1]


def test_mad_outlier_with_zero_dispersion_flags_nothing():
    """A constant vector has MAD 0; the rule must not divide by zero."""
    assert not mad_outlier(np.full(50, 7.0), n_mads=3.0).any()


def test_zebrafish_mito_prefix_is_lowercase():
    """The human 'MT-' convention matches no zebrafish gene and fails silently."""
    adata = _toy_adata()
    annotate_gene_classes(adata, QCConfig(mito_prefix="mt-"))
    assert adata.var["mt"].sum() == 2

    adata_human_prefix = _toy_adata()
    annotate_gene_classes(adata_human_prefix, QCConfig(mito_prefix="zz-"))
    assert adata_human_prefix.var["mt"].sum() == 0
    assert (adata_human_prefix.obs["pct_counts_mt"] == 0).all()


def test_flag_low_quality_catches_high_mito_cell():
    adata = _toy_adata()
    adata.X[0, :2] = 50_000.0  # cell 0 becomes almost entirely mitochondrial
    annotate_gene_classes(adata, QCConfig())
    flag_low_quality(adata, QCConfig(min_genes_per_cell=0))
    assert bool(adata.obs["qc_fail"].iloc[0])


def test_flag_low_quality_is_not_indiscriminate():
    adata = _toy_adata()
    annotate_gene_classes(adata, QCConfig())
    flag_low_quality(adata, QCConfig(min_genes_per_cell=0))
    assert adata.obs["qc_fail"].mean() < 0.2


def _multi_batch_adata(n_per_batch: int = 400, n_batches: int = 3,
                       seed: int = 0) -> ad.AnnData:
    """Several batches, so Scrublet runs its per-batch path as it does on real data."""
    rng = np.random.default_rng(seed)
    genes = ["mt-nd1", "mt-co1"] + [f"gene{i:03d}" for i in range(298)]
    blocks, samples = [], []
    for batch in range(n_batches):
        blocks.append(rng.negative_binomial(8, 0.3,
                                            size=(n_per_batch, len(genes))).astype(float))
        samples += [f"lib{batch}"] * n_per_batch
    adata = ad.AnnData(np.vstack(blocks))
    adata.var_names = genes
    adata.obs["sample"] = samples
    adata.obs_names = [f"{s}_{i:05d}" for i, s in enumerate(samples)]
    return adata


def test_doublet_detection_survives_barcodes_stripped_to_zero():
    """Regression: CellBender can leave a barcode with no counts at all.

    Scrublet's own preprocessing drops cells below three counts, after which
    scanpy cannot reindex the per-batch results onto the full set and raises
    `KeyError: [...] not in index`. Before the fix that exception was caught and
    the run continued with every `predicted_doublet` False — producing a
    complete-looking result set from which no doublet had been removed. Cell
    Ranger never emits such barcodes, so only ambient-corrected input hits this.
    """
    adata = _multi_batch_adata()
    # Strip a handful of barcodes across different batches, exactly as ambient
    # correction does: zero counts, and one to two counts either side of
    # Scrublet's internal threshold.
    stripped = [0, 1, 400, 401, 800]
    adata.X[stripped, :] = 0.0
    adata.X[401, 0] = 2.0

    annotate_gene_classes(adata, QCConfig())
    flag_low_quality(adata, QCConfig())
    status = detect_doublets(adata, QCConfig(), seed=0)

    assert status.startswith("ok:"), f"doublet detection did not run: {status}"
    # The stripped barcodes fail QC on gene count, so they must be excluded from
    # scoring rather than fed to Scrublet.
    assert adata.obs["qc_fail"].to_numpy()[stripped].all()
    assert np.isnan(adata.obs["doublet_score"].to_numpy()[stripped]).all()
    # QC-passing cells were genuinely scored.
    passing = ~adata.obs["qc_fail"].to_numpy()
    assert not np.isnan(adata.obs["doublet_score"].to_numpy()[passing]).all()


def test_degenerate_automatic_threshold_is_refused():
    """Regression: a threshold no simulated doublet reaches must not pass.

    On GSE158142 Scrublet returned 0.825 for both timepoints — above the highest
    observed cell and above 100% of its own simulated doublets. Zero cells were
    called, in every run, and the result was indistinguishable from a dataset
    with no doublets. Reported as a finding it would have been simply false.
    """
    adata = _multi_batch_adata(n_per_batch=60, n_batches=2)
    adata.obs["doublet_score"] = np.linspace(0.0, 0.7, adata.n_obs)
    adata.obs["predicted_doublet"] = False
    # Mirrors the real shape: a few simulated doublets DO exceed the threshold
    # (12 of 102,462 on GSE158142), so a check that only asks whether the
    # simulated set clears it passes while no real cell ever can.
    sim = np.concatenate([np.linspace(0.002, 0.72, 500), [0.83, 0.84, 0.91]])
    adata.uns["scrublet"] = {"batches": {
        "lib0": {"threshold": 0.825, "doublet_scores_sim": sim},
        "lib1": {"threshold": 0.828, "doublet_scores_sim": sim},
    }}
    with pytest.raises(ValueError, match="degenerate"):
        _apply_threshold(adata, QCConfig(doublet_threshold="auto"))


def test_automatic_threshold_accepted_when_it_can_call_something():
    adata = _multi_batch_adata(n_per_batch=60, n_batches=1)
    adata.obs["doublet_score"] = np.linspace(0.0, 0.9, adata.n_obs)
    adata.uns["scrublet"] = {
        "threshold": 0.4,
        "doublet_scores_sim": np.linspace(0.002, 0.95, 500),
    }
    note = _apply_threshold(adata, QCConfig(doublet_threshold="auto"))
    assert "auto" in note


def test_explicit_threshold_overrides_and_reports_sensitivity():
    """Sensitivity must be reported: the call rate alone reads as a doublet rate."""
    adata = _multi_batch_adata(n_per_batch=100, n_batches=1)
    adata.obs["doublet_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.obs["predicted_doublet"] = False
    adata.uns["scrublet"] = {
        "threshold": 0.825,  # the degenerate one; the explicit value must win
        "doublet_scores_sim": np.linspace(0.0, 1.0, 1000),
    }
    note = _apply_threshold(adata, QCConfig(doublet_threshold=0.25))

    called = adata.obs["predicted_doublet"].to_numpy()
    assert called.sum() > 0, "an explicit threshold below the score range must call cells"
    assert (adata.obs["doublet_score"].to_numpy()[called] > 0.25).all()
    assert "explicit" in note and "sensitivity" in note


def test_doublet_detection_status_distinguishes_disabled_from_found_none():
    """A run that skipped the step must not look like one that found nothing."""
    adata = _multi_batch_adata(n_per_batch=50, n_batches=2)
    annotate_gene_classes(adata, QCConfig())
    flag_low_quality(adata, QCConfig())
    status = detect_doublets(adata, QCConfig(run_doublet_detection=False), seed=0)
    assert status == "disabled"
    assert not adata.obs["predicted_doublet"].any()


def test_config_reports_missing_input_files(tmp_path):
    cfg = Config(
        samples=[SampleSpec(name="day5", h5=tmp_path / "absent.h5", timepoint="5dpf")],
        outdir=tmp_path / "out",
        markers_file=tmp_path / "markers_ocular.yaml",
        compartment_markers=tmp_path / "markers_compartment.yaml",
    )
    with pytest.raises(FileNotFoundError, match="absent.h5"):
        cfg.validate_inputs()


def test_resolve_makes_relative_paths_absolute(tmp_path):
    """A notebook and the CLI have different working directories; one config
    must serve both without carrying machine-specific absolute paths."""
    from pathlib import Path

    cfg = Config(
        samples=[SampleSpec(name="day5", h5=Path("data/day5/x.h5"), timepoint="5dpf"),
                 SampleSpec(name="day8", h5=Path("data/day8/x.h5"), timepoint="8dpf")],
        outdir=Path("results"),
        markers_file=Path("config/markers_ocular.yaml"),
        compartment_markers=Path("config/markers_compartment.yaml"),
    )
    resolved = cfg.resolve(tmp_path)

    assert all(s.h5.is_absolute() for s in resolved.samples)
    assert resolved.samples[0].h5 == tmp_path.resolve() / "data/day5/x.h5"
    assert resolved.outdir == tmp_path.resolve() / "results"
    assert resolved.markers_file == tmp_path.resolve() / "config/markers_ocular.yaml"
    assert resolved.compartment_markers == tmp_path.resolve() / "config/markers_compartment.yaml"
    # Non-destructive: the original is untouched.
    assert not cfg.outdir.is_absolute()


def test_resolve_leaves_absolute_paths_alone(tmp_path):
    from pathlib import Path

    absolute = tmp_path / "elsewhere" / "day5.h5"
    cfg = Config(
        samples=[SampleSpec(name="day5", h5=absolute, timepoint="5dpf")],
        outdir=Path("results"),
        markers_file=Path("config/markers_ocular.yaml"),
        compartment_markers=Path("config/markers_compartment.yaml"),
    )
    assert cfg.resolve("/some/other/root").samples[0].h5 == absolute


def test_cli_dry_run_validates_without_running(tmp_path, capsys):
    """`--dry-run` must confirm the plan and produce no outputs."""
    import subprocess
    import sys as _sys

    from screye.cli import main

    root = Path(__file__).resolve().parents[1]
    data = tmp_path / "data"
    subprocess.run(
        [_sys.executable, str(root / "tests/fixtures/make_synthetic_h5.py"),
         "--markers", str(root / "config/markers_compartment.yaml"),
         "--outdir", str(data), "--n-cells", "50"],
        check=True, capture_output=True,
    )

    cfg_text = (root / "config/config.yaml").read_text()
    cfg_text = cfg_text.replace("data/day5_filtered_feature_bc_matrix.h5",
                                str(data / "day5_filtered_feature_bc_matrix.h5"))
    cfg_text = cfg_text.replace("data/day8_filtered_feature_bc_matrix.h5",
                                str(data / "day8_filtered_feature_bc_matrix.h5"))
    cfg_text = cfg_text.replace("markers_file: config/markers_ocular.yaml",
                                f"markers_file: {root / 'config/markers_ocular.yaml'}")
    cfg_text = cfg_text.replace("compartment_markers: config/markers_compartment.yaml",
                                f"compartment_markers: {root / 'config/markers_compartment.yaml'}")
    cfg_text = cfg_text.replace("outdir: results", f"outdir: {tmp_path / 'results'}")
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(cfg_text)

    assert main(["-c", str(cfg_path), "--root", str(root), "--dry-run"]) == 0
    assert "day5" in capsys.readouterr().out
    assert not (tmp_path / "results").exists()


def test_cli_reports_missing_config_without_traceback(tmp_path):
    from screye.cli import main

    assert main(["-c", str(tmp_path / "absent.yaml"),
                 "--root", str(tmp_path), "--dry-run"]) == 1
