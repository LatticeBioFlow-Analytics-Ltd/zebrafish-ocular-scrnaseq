"""Unit tests for two-pass compartment assignment and publication figure sizing.

These cover the decisions that are easy to get wrong quietly: an assignment made
on a margin too narrow to justify it, a subset that silently inherits stale
clustering state from pass 1, and a figure that claims to be journal-sized but
is not.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from screye import figures as figs
from screye.compartments import assign_by_margin, load_panel

REPO = Path(__file__).resolve().parents[1]


# --- compartment assignment -------------------------------------------------

def test_clear_winner_is_assigned_without_a_flag():
    scores = pd.DataFrame(
        {"ocular_photoreceptor": [2.0], "neuron_differentiated": [0.1]},
        index=["0"],
    )
    labels, margins = assign_by_margin(scores, min_margin=0.05)
    assert labels["0"] == "ocular_photoreceptor"
    assert margins["0"] == pytest.approx(1.9)
    assert "low confidence" not in labels["0"]


def test_narrow_margin_is_flagged_rather_than_asserted():
    """Two near-equal scores mean the panel cannot separate those identities.
    Reporting the winner silently would be an unearned claim."""
    scores = pd.DataFrame(
        {"ocular_retinal_neuron_glia": [1.00], "neuron_differentiated": [0.99]},
        index=["7"],
    )
    labels, _ = assign_by_margin(scores, min_margin=0.05)
    assert "low confidence" in labels["7"]
    assert labels["7"].startswith("ocular_retinal_neuron_glia")


def test_margin_threshold_is_respected_at_the_boundary():
    scores = pd.DataFrame({"a": [1.0], "b": [0.94]}, index=["0"])
    assert "low confidence" not in assign_by_margin(scores, 0.05)[0]["0"]
    assert "low confidence" in assign_by_margin(scores, 0.10)[0]["0"]


# --- marker panels ----------------------------------------------------------

def test_both_shipped_panels_parse_and_are_non_empty():
    for name in ("markers_compartment.yaml", "markers_ocular.yaml"):
        panel = load_panel(REPO / "config" / name)
        assert panel, f"{name} is empty"
        assert all(genes for genes in panel.values()), f"{name} has an empty set"


def test_compartment_panel_covers_the_ocular_prefix():
    """The coarse grouping in assign_compartments keys off an 'ocular_' prefix,
    so the panel and that logic must not drift apart."""
    panel = load_panel(REPO / "config" / "markers_compartment.yaml")
    assert any(k.startswith("ocular_") for k in panel)


def test_ocular_panel_resolves_retinal_classes_the_broad_panel_does_not():
    """The point of pass 2 is finer resolution; if the ocular panel were not
    more specific than the compartment panel, the second pass would be
    pointless."""
    broad = load_panel(REPO / "config" / "markers_compartment.yaml")
    fine = load_panel(REPO / "config" / "markers_ocular.yaml")
    for cell_type in ("rod_photoreceptor", "cone_photoreceptor", "bipolar_cell",
                      "amacrine_cell", "horizontal_cell", "retinal_ganglion",
                      "muller_glia"):
        assert cell_type in fine, f"{cell_type} missing from the ocular panel"
        assert cell_type not in broad


# --- publication figures ----------------------------------------------------

def test_diverging_colormap_puts_red_at_the_top_and_blue_at_the_bottom():
    """The project-wide convention: increased = red, decreased = blue. If this
    ever inverts, every figure silently reverses meaning."""
    low = figs.DIVERGING(0.0)
    high = figs.DIVERGING(1.0)
    assert low[2] > low[0], "low end of the scale should be blue"
    assert high[0] > high[2], "high end of the scale should be red"


def test_diverging_norm_centres_white_on_zero():
    """With an asymmetric range, an uncentred norm would render genuinely
    unchanged features as pale blue - a visual claim of decrease."""
    norm = figs.diverging_norm(-1.0, 4.0)
    assert norm(0.0) == pytest.approx(0.5)
    assert norm(-1.0) == pytest.approx(0.0)
    assert norm(4.0) == pytest.approx(1.0)


def test_fit_to_width_converges_on_the_journal_width(tmp_path):
    """scanpy grids arrive several hundred mm wide; the fit must bring them to
    the requested width despite legends that sit outside the axes."""
    figs.apply_style()
    fig, ax = plt.subplots(figsize=(15, 8))
    ax.scatter(np.random.default_rng(0).normal(size=50),
               np.random.default_rng(1).normal(size=50), label="a long legend entry")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))

    achieved = figs.fit_to_width(fig, figs.DOUBLE_COL)
    assert achieved == pytest.approx(figs.DOUBLE_COL, rel=0.02)
    plt.close(fig)


def test_save_figure_writes_vector_formats_by_default(tmp_path):
    """DMM accepts EPS, PDF and SVG for line art and warns that JPEG/TIFF may
    delay production, so a vector file must always be produced - a PNG-only
    run would leave nothing submittable."""
    figs.apply_style()
    plt.subplots()
    written = figs.save_figure(tmp_path, "test_fig")
    suffixes = {p.suffix.lstrip(".") for p in written}
    assert suffixes >= set(figs.VECTOR_FORMATS)
    assert all(p.exists() and p.stat().st_size > 0 for p in written)


def test_vector_output_keeps_text_as_text(tmp_path):
    """DMM: save text as 'text', not 'curves' or 'outlines'. If labels are
    converted to paths the journal cannot edit them during production."""
    figs.apply_style()
    _fig, ax = plt.subplots()
    ax.set_xlabel("a distinctive axis label")
    figs.save_figure(tmp_path, "vec", formats=("svg",))
    assert "a distinctive axis label" in (tmp_path / "vec.svg").read_text()


def test_type_sizes_match_the_journal_specification():
    """8 pt for labelling, 12 pt bold for panel letters."""
    import matplotlib as mpl

    figs.apply_style()
    assert figs.LABEL_PT == 8
    assert figs.PANEL_LETTER_PT == 12
    for key in ("font.size", "axes.labelsize", "xtick.labelsize",
                "ytick.labelsize", "legend.fontsize"):
        assert mpl.rcParams[key] == figs.LABEL_PT, key


def test_panel_letters_are_upper_case_and_bold():
    figs.apply_style()
    fig, ax = plt.subplots()
    figs.panel_label(ax, "a")
    text = ax.texts[-1]
    assert text.get_text() == "A"
    assert text.get_fontweight() == "bold"
    assert text.get_fontsize() == figs.PANEL_LETTER_PT
    plt.close(fig)


def test_page_limits_match_the_journal_specification():
    """180 mm x 210 mm including lettering and labels."""
    assert figs.MAX_WIDTH / figs.MM == pytest.approx(180)
    assert figs.MAX_HEIGHT / figs.MM == pytest.approx(210)
    assert figs.DOUBLE_COL <= figs.MAX_WIDTH


def test_working_height_leaves_room_for_a_legend():
    """DMM's page schematics set the legend in its own box beneath the figure,
    so a figure filling the full 210 mm leaves nothing for it."""
    assert figs.WORKING_HEIGHT < figs.MAX_HEIGHT


def test_intermediate_column_widths_are_flagged(caplog):
    """DMM prefers single or double column and shows a 1.5-column layout as a
    negative example."""
    with caplog.at_level("WARNING"):
        figs.assert_column_width(120 * figs.MM, "odd")
    assert "neither a single" in caplog.text

    caplog.clear()
    with caplog.at_level("WARNING"):
        figs.assert_column_width(figs.SINGLE_COL, "single")
        figs.assert_column_width(figs.DOUBLE_COL, "double")
    assert caplog.text == ""


def test_background_fraction_detects_sparse_layouts():
    """A proxy for DMM's "maximise data, minimise background". A single small
    axes in a large canvas should score far higher than one that fills it."""
    figs.apply_style()
    sparse, ax = plt.subplots(figsize=(6, 6))
    ax.set_position([0.1, 0.1, 0.2, 0.2])
    full, ax2 = plt.subplots(figsize=(6, 6))
    ax2.set_position([0.0, 0.0, 1.0, 1.0])

    assert figs.background_fraction(sparse) > figs.background_fraction(full)
    assert figs.background_fraction(full) < 0.1
    plt.close("all")


def test_saved_figure_is_within_the_double_column_limit(tmp_path):
    """The end-to-end guarantee: what lands on disk is the width claimed."""
    from PIL import Image

    figs.apply_style()
    _fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot([0, 1], [0, 1], label="series")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    figs.save_figure(tmp_path, "wide", formats=("png",), dpi=300)

    with Image.open(tmp_path / "wide.png") as img:   # context-managed: an open
        width_mm = img.size[0] / 300 * 25.4          # handle trips the strict
                                                     # warning filter in CI
    assert width_mm <= figs.MAX_WIDTH / figs.MM, (
        f"saved figure is {width_mm:.1f} mm wide, over the 180 mm limit")


# --- panel spacing ----------------------------------------------------------

def _obs_frame(**columns):
    import anndata as ad_
    import numpy as np_
    import pandas as pd_

    n = len(next(iter(columns.values())))
    a = ad_.AnnData(np_.zeros((n, 2), dtype="float32"))
    for key, values in columns.items():
        a.obs[key] = pd_.Categorical(values)
    return a


def test_spacing_widens_for_many_legend_entries():
    """scanpy draws each legend in the gap right of its panel, so a 16-level
    cluster key needs more room than a two-level sample key."""
    few = _obs_frame(sample=["day5", "day8"] * 50)
    many = _obs_frame(leiden=[str(i % 16) for i in range(100)])
    assert figs.panel_spacing_for(many, ["leiden"]) > figs.panel_spacing_for(few, ["sample"])


def test_spacing_widens_for_long_labels():
    short = _obs_frame(k=["a", "b"] * 50)
    long = _obs_frame(k=["cmz_progenitor (low confidence)", "photoreceptor_precursor"] * 50)
    assert figs.panel_spacing_for(long, ["k"]) > figs.panel_spacing_for(short, ["k"])


def test_spacing_is_bounded():
    """Legends must not be allowed to consume the whole page."""
    extreme = _obs_frame(k=[f"a_very_long_cell_type_label_{i}" for i in range(200)])
    spacing = figs.panel_spacing_for(extreme, ["k"])
    assert 0.30 <= spacing <= 0.85


def test_spacing_ignores_missing_and_numeric_keys():
    a = _obs_frame(sample=["day5", "day8"] * 50)
    a.obs["counts"] = list(range(100))
    assert figs.panel_spacing_for(a, ["absent", "counts"]) == pytest.approx(0.30)


def test_dot_sizing_is_reduced_from_the_scanpy_default():
    """scanpy's 200 pt^2 default leaves dots touching at 180 mm width."""
    from scanpy.plotting import DotPlot

    from screye import pipeline

    assert pipeline.LARGEST_DOT < DotPlot.DEFAULT_LARGEST_DOT
