"""Execute the walkthrough notebook end to end.

A notebook that does not run is worse than no notebook: it is the first thing a
reviewer opens and the first thing that can embarrass you. This test catches
imports dropped during editing, hardcoded relative paths, and cells left in an
order that only worked in the author's kernel.

Marked `integration` because it runs the full analysis (~2 min).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
NOTEBOOK = REPO / "notebooks" / "01_walkthrough.ipynb"

pytestmark = pytest.mark.integration


def test_walkthrough_notebook_runs(tmp_path, monkeypatch):
    nbformat = pytest.importorskip("nbformat")
    nbclient = pytest.importorskip("nbclient")
    yaml = pytest.importorskip("yaml")

    # Point the config at synthetic data so the test needs no real libraries.
    # The notebook reads SCREYE_CONFIG: without it, it would silently analyse
    # whatever real data sits in data/ - unbounded runtime and a test that
    # passes or times out depending on which machine it runs on.
    data = tmp_path / "data"
    subprocess.run(
        [sys.executable, str(REPO / "tests/fixtures/make_synthetic_h5.py"),
         "--markers", str(REPO / "config/markers_compartment.yaml"),
         "--outdir", str(data), "--n-cells", "300"],
        check=True, capture_output=True,
    )
    base = yaml.safe_load((REPO / "config/config.yaml").read_text())
    base["samples"] = [
        {"name": "day5", "h5": str(data / "day5_filtered_feature_bc_matrix.h5"),
         "timepoint": "5dpf"},
        {"name": "day8", "h5": str(data / "day8_filtered_feature_bc_matrix.h5"),
         "timepoint": "8dpf"},
    ]
    base["outdir"] = str(tmp_path / "results")
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(base))

    monkeypatch.setenv("PYTHONPATH", str(REPO / "src"))
    monkeypatch.setenv("SCREYE_CONFIG", str(cfg_path))

    nb = nbformat.read(NOTEBOOK, as_version=4)
    client = nbclient.NotebookClient(
        nb, timeout=900, kernel_name="python3",
        resources={"metadata": {"path": str(NOTEBOOK.parent)}},
    )
    client.execute()  # raises CellExecutionError on any failing cell


def test_notebook_contains_no_hardcoded_relative_paths():
    """Paths must come from Config.resolve(), not from '../' guesses."""
    import json

    nb = json.loads(NOTEBOOK.read_text())
    offenders = [
        i for i, cell in enumerate(nb["cells"])
        if cell["cell_type"] == "code" and "../" in "".join(cell["source"])
    ]
    assert not offenders, f"hardcoded relative paths in cells {offenders}"
