"""Unit tests for the article-panel detection statistics.

The figures are exercised by the integration run; what needs pinning here is
the arithmetic the acceptance criteria depend on - detection fractions per
cluster and stage - because an off-by-one in a boolean mask produces plausible
wrong percentages that no rendering check would catch.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from screye.article import _acceptance_criteria, _rpe_mask, detection_table

SPEC = {
    "rpe_cell_type": "rpe",
    "rpe_identity_markers": ["ident"],
    "panel": {
        "RPE identity": ["ident"],
        "abundance-matched positive": ["matched"],
        "target": ["target1"],
        "ambient sentinel": ["rho"],
        "out-of-tissue": ["absent_gene"],
    },
}


@pytest.fixture
def toy() -> ad.AnnData:
    """20 cells: 10 RPE (first half 5dpf), 10 other; hand-set detection.

    target1: detected in 8/10 RPE cells and 1/10 others.
    rho:     detected in 2/10 RPE cells and 9/10 others.
    """
    n = 20
    genes = ["ident", "matched", "target1", "rho"]
    x = np.zeros((n, len(genes)), dtype=np.float32)
    x[:10, genes.index("target1")][:8] = 1.0
    x[10:, genes.index("target1")][:1] = 1.0
    x[:10, genes.index("rho")][:2] = 1.0
    x[10:, genes.index("rho")][:9] = 1.0
    x[:, genes.index("ident")] = 1.0
    x[:10, genes.index("matched")][:7] = 1.0
    a = ad.AnnData(x)
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(["rpe"] * 8 + ["rpe (low confidence)"] * 2
                                        + ["other"] * 10)
    a.obs["timepoint"] = pd.Categorical(["5dpf"] * 5 + ["8dpf"] * 5
                                        + ["5dpf"] * 5 + ["8dpf"] * 5)
    a.raw = a
    return a


def test_rpe_mask_includes_low_confidence(toy):
    assert _rpe_mask(toy, "rpe").sum() == 10


def test_detection_fractions_match_hand_counts(toy):
    table = detection_table(toy, SPEC)
    pooled = table[table["stage"] == "pooled"].set_index("gene")
    assert pooled.loc["target1", "pct_rpe"] == pytest.approx(80.0)
    assert pooled.loc["target1", "pct_non_rpe"] == pytest.approx(10.0)
    assert pooled.loc["rho", "pct_rpe"] == pytest.approx(20.0)
    # per-stage split: all 8 target1-positive RPE cells were laid out first,
    # so 5dpf has 5/5 and 8dpf 3/5
    d5 = table[table["stage"] == "5dpf"].set_index("gene")
    d8 = table[table["stage"] == "8dpf"].set_index("gene")
    assert d5.loc["target1", "pct_rpe"] == pytest.approx(100.0)
    assert d8.loc["target1", "pct_rpe"] == pytest.approx(60.0)


def test_missing_gene_is_reported_not_silently_dropped(toy):
    table = detection_table(toy, SPEC)
    row = table[table["gene"] == "absent_gene"]
    assert len(row) == 1
    assert not row["detected_in_dataset"].iloc[0]
    assert np.isnan(row["pct_rpe"].iloc[0])


def test_acceptance_criteria_verdicts(toy):
    criteria = _acceptance_criteria(detection_table(toy, SPEC), SPEC)
    by_number = {int(c.split(".")[0]): v
                 for c, v in zip(criteria["criterion"], criteria["verdict"])}
    assert by_number[1] == "pass"          # target 80% > rho 20% in RPE
    assert by_number[3] == "pass"          # 80% RPE vs 10% elsewhere
    assert by_number[4] == "pending"       # ambient correction is external
    assert by_number[5].startswith("FAIL")  # 5 RPE cells per stage < 50
