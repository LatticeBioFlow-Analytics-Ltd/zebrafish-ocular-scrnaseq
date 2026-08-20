"""Sample-sheet expansion: the handover from the scRNA_analysis workflow.

The upstream Snakemake project publishes one ambient-corrected h5 per 10x
library and keeps the library-to-timepoint mapping in config/libraries.tsv.
These tests pin the contract: the sheet is the source of truth, nothing is
parsed from filenames, and misdeclared configs fail loudly at load time.
"""

from __future__ import annotations

import pytest
import yaml

from screye.config import Config, SampleSheet

SHEET = (
    "library\ttimepoint\tn_runs\texpect_cells\tdescription\n"
    "zBr5dpf1_S7\t5dpf\t4\t\t\n"
    "zBr5dpf2_S3\t5dpf\t4\t\t\n"
    "zBr8dpf1_S2\t8dpf\t8\t\t\n"
)


def _write_yaml(tmp_path, body: dict) -> str:
    base = {
        "outdir": "results",
        "markers_file": "config/markers_ocular.yaml",
        "compartment_markers": "config/markers_compartment.yaml",
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({**base, **body}))
    return path


def test_sheet_expands_to_one_spec_per_library(tmp_path):
    (tmp_path / "sheet.tsv").write_text(SHEET)
    cfg_path = _write_yaml(tmp_path, {
        "sample_sheet": {
            "libraries": "sheet.tsv",
            "h5_pattern": "h5/{library}_cellbender_filtered.h5",
        },
    })

    cfg = Config.from_yaml(cfg_path).resolve(tmp_path)

    assert [s.name for s in cfg.samples] == [
        "zBr5dpf1_S7", "zBr5dpf2_S3", "zBr8dpf1_S2"]
    assert [s.timepoint for s in cfg.samples] == ["5dpf", "5dpf", "8dpf"]
    # h5 paths come from the pattern, resolved against the root - not parsed
    # from anything.
    assert cfg.samples[0].h5 == tmp_path / "h5/zBr5dpf1_S7_cellbender_filtered.h5"
    assert all(s.h5.is_absolute() for s in cfg.samples)


def test_sheet_and_samples_together_is_an_error(tmp_path):
    (tmp_path / "sheet.tsv").write_text(SHEET)
    cfg_path = _write_yaml(tmp_path, {
        "samples": [{"name": "day5", "h5": "a.h5", "timepoint": "5dpf"}],
        "sample_sheet": {"libraries": "sheet.tsv", "h5_pattern": "{library}.h5"},
    })
    with pytest.raises(ValueError, match="exactly one"):
        Config.from_yaml(cfg_path)


def test_neither_samples_nor_sheet_is_an_error(tmp_path):
    cfg_path = _write_yaml(tmp_path, {})
    with pytest.raises(ValueError, match="exactly one"):
        Config.from_yaml(cfg_path)


def test_pattern_without_placeholder_is_rejected():
    with pytest.raises(ValueError, match="placeholder"):
        SampleSheet(libraries="sheet.tsv", h5_pattern="fixed_name.h5")


def test_sheet_missing_timepoint_column_is_rejected(tmp_path):
    (tmp_path / "sheet.tsv").write_text("library\tstage\nzBr5dpf1_S7\t5dpf\n")
    cfg_path = _write_yaml(tmp_path, {
        "sample_sheet": {"libraries": "sheet.tsv", "h5_pattern": "{library}.h5"},
    })
    with pytest.raises(ValueError, match="timepoint"):
        Config.from_yaml(cfg_path).resolve(tmp_path)


def test_unexpanded_sheet_fails_validation_with_guidance(tmp_path):
    (tmp_path / "sheet.tsv").write_text(SHEET)
    cfg_path = _write_yaml(tmp_path, {
        "sample_sheet": {"libraries": "sheet.tsv", "h5_pattern": "{library}.h5"},
    })
    cfg = Config.from_yaml(cfg_path)  # .resolve() deliberately not called
    with pytest.raises(ValueError, match="resolve"):
        cfg.validate_inputs()


def test_real_snakemake_sheet_expands_if_present():
    """Against the actual sibling scRNA_analysis project, when it exists.

    Guards the real contract (column names, 22 libraries, 9+13 split) rather
    than a fixture's idea of it. Skipped cleanly on machines without the
    upstream checkout.
    """
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    sheet = repo.parent / "scRNA_analysis/config/libraries.tsv"
    if not sheet.exists():
        pytest.skip("scRNA_analysis project not checked out beside this repo")

    specs = SampleSheet(
        libraries=sheet,
        h5_pattern="../scRNA_analysis/results/h5/cellbender/{library}_cellbender_filtered.h5",
    ).expand(repo)

    timepoints = [s.timepoint for s in specs]
    assert len(specs) == 22
    assert timepoints.count("5dpf") == 9
    assert timepoints.count("8dpf") == 13
