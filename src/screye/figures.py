"""Publication-ready figure styling for Disease Models & Mechanisms.

Encodes the figure requirements at
https://journals.biologists.com/dmm/pages/manuscript-prep (section 4, last
updated 22 July 2026). Re-check before submission; journal specifications
change and this file is a convenience, not an authority.

Specification as implemented
----------------------------
  maximum size      180 mm x 210 mm including lettering and labels
  resolution        300 ppi photographic; 600 ppi where photographs and line
                    art or text are combined
  panel letters     12 pt bold uppercase (A, B, C)
  other labelling   8 pt Arial, sentence case; headings bold
  fonts             Arial or Helvetica throughout, saved as text not outlines
  colour            RGB, no spot/Pantone colours, no embedded colour profile
  line art formats  EPS, PDF, SVG. JPEG and TIFF for line art may delay
                    production, so PNG here is for review and preprints only.

Colour convention
-----------------
Diverging quantities (log fold change, z-scored expression, enrichment) use RED
for positive/increased and BLUE for negative/decreased, centred on zero, applied
consistently across every figure in the project.

This is compatible with, and motivated by, DMM's colour-vision guidance: the
journal asks authors to avoid red/green for two-channel displays and encourages
colour-blind-safe choices particularly in heatmaps. Red-to-blue is the
recommended alternative axis rather than the prohibited one. The palette is
ColorBrewer RdBu reversed. Sequential quantities use viridis (perceptually
uniform, greyscale-safe); categorical variables use Okabe-Ito.

Not verified
------------
DMM links a separate figure-layout PDF from section 4.1 which could not be
retrieved when this module was written. It may carry panel-arrangement guidance
not reflected here.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

log = logging.getLogger(__name__)

# --- dimensions (mm -> inches) ---------------------------------------------
MM = 1 / 25.4

# DMM: "The maximum figure size, including lettering and labels, is
# 180 mm x 210 mm." Both are hard ceilings, not targets - the journal also asks
# for "the smallest size that will convey the essential scientific information".
MAX_WIDTH = 180 * MM
MAX_HEIGHT = 210 * MM
# Fit target sits just inside the ceiling. Rasterising at 600 dpi rounds to
# whole pixels, and savefig recomputes the tight bounding box after the fit, so
# aiming at exactly 180 mm lands a fraction over. 0.5 mm of headroom is
# invisible in print and keeps every output provably inside the limit.
SAFETY_MARGIN = 0.5 * MM
DOUBLE_COL = MAX_WIDTH - SAFETY_MARGIN

# DMM's image guidelines show the journal page as two text columns into which
# figures fit as either one column ("single column / half-page width") or both.
# Half of 180 mm less a gutter gives 85 mm. Intermediate widths are shown as a
# negative example: a 1.5-column layout is criticised for suboptimal layout and
# unminimised background.
SINGLE_COL = 85 * MM
COLUMN_WIDTHS = (SINGLE_COL, DOUBLE_COL)

# In the page schematics the legend sits in its own box beneath the figure, so
# a figure that fills the full 210 mm leaves no room for its own legend. 210 mm
# remains the hard ceiling; this is the height to design to.
WORKING_HEIGHT = 175 * MM

# 600 ppi is DMM's preference for figures combining photographs with line art or
# text, and is a safe default for everything raster. Vector output has no dpi.
DPI_SUBMISSION = 600
DPI_REVIEW = 300

# Type sizes, from DMM section 4.3.1.
PANEL_LETTER_PT = 12   # "12 pt bold uppercase letters (A, B, C, etc.)"
LABEL_PT = 8           # "Other labelling should be 8 pt Arial font"

# DMM accepts EPS, PDF and SVG for line art and warns that JPEG/TIFF "may delay
# production". PNG is included for review, preprints and README display only.
VECTOR_FORMATS = ("pdf", "svg")
DEFAULT_FORMATS = ("png", "pdf", "svg")
COMPLIANT_FONTS = ("Arial", "Helvetica")

# --- colour convention ------------------------------------------------------
# ColorBrewer RdBu, reversed so that high = red. Defined explicitly rather than
# via plt.get_cmap("RdBu_r") so the exact hex values are recorded in the
# repository and cannot shift with a matplotlib version change.
_RDBU = [
    "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
    "#f7f7f7",
    "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f",
]
DIVERGING = LinearSegmentedColormap.from_list("screye_rdbu", _RDBU, N=256)
SEQUENTIAL = "viridis"

# Categorical palette for cell types: Okabe-Ito, designed for colour-vision
# deficiency, extended with greys for overflow categories.
CATEGORICAL = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#000000", "#999999", "#7F3C8D",
    "#11A579", "#3969AC", "#F2B701", "#E73F74", "#80BA5A",
    "#E68310", "#008695", "#CF1C90", "#f97b72", "#4b4b8f",
]


# --- content-driven sizing --------------------------------------------------
# Grid figures (dotplots, heatmaps) are sized from their content rather than
# drawn at scanpy's default and rescaled afterwards. Fonts and scatter markers
# are specified in points, so rescaling the canvas changes their size relative
# to everything positioned in axes coordinates: shrinking a 400 mm dotplot to
# 180 mm leaves the 8 pt labels and 45 pt^2 dots at full size while the space
# between them halves, which is exactly the overlapping-dots, colliding-legends
# failure mode. Choosing the canvas from the number of rows and columns at the
# final printed size means nothing needs rescaling and nothing collides.
#
# Pitch: an 8 pt label is ~2.8 mm tall, so 3.3 mm per row/column keeps adjacent
# rotated gene labels apart and leaves the largest dot (sqrt(45) ~ 6.7 pt
# ~ 2.4 mm diameter) clear of its neighbours.
DOT_PITCH = 0.13          # inches per dotplot row or column
# Width scanpy reserves for the colourbar + size legend, in inches. 1.6 rather
# than the default 1.5: the colourbar title and the size key collide otherwise.
DOT_LEGEND_WIDTH = 1.6
# Allowance for the dendrogram panel and row-label overhang beside the grid.
DOT_EXTRAS_ALLOWANCE = 1.1
# Columns that fit a double-column figure at DOT_PITCH once the legend block
# and the dendrogram have taken their share. Chunk any wider plot.
MAX_DOT_COLUMNS = int((DOUBLE_COL - DOT_LEGEND_WIDTH - DOT_EXTRAS_ALLOWANCE)
                      / DOT_PITCH)


def timepoint_shades(sample_to_timepoint: dict[str, str]) -> dict[str, tuple]:
    """One colour per library, shaded light-to-dark within its timepoint.

    With one library per timepoint this reduces to the first Okabe-Ito colours.
    With many libraries a flat categorical palette makes the stacked composition
    bars unreadable (22 unrelated hues) and hides the structure that matters:
    which timepoint a library belongs to. Shades of one hue per timepoint keep
    both levels legible.
    """
    from matplotlib.colors import to_rgb

    order = list(dict.fromkeys(sample_to_timepoint.values()))
    shades: dict[str, tuple] = {}
    for j, tp in enumerate(order):
        libs = [s for s, t in sample_to_timepoint.items() if t == tp]
        base = np.asarray(to_rgb(CATEGORICAL[j % len(CATEGORICAL)]))
        for k, lib in enumerate(libs):
            # Blend toward white; never lighter than 45% strength so the
            # palest shade still reads against a white background.
            w = 1.0 if len(libs) == 1 else 0.45 + 0.55 * (k + 1) / len(libs)
            shades[lib] = tuple(1.0 - w * (1.0 - base))
    return shades


def diverging_norm(vmin: float, vmax: float) -> TwoSlopeNorm:
    """Normalisation anchored at zero, so white always means "no change".

    Without this, a range of e.g. -1 to +4 puts white at +1.5 and every genuinely
    unchanged feature is rendered pale blue - a visual claim of decrease where
    the data show none.
    """
    return TwoSlopeNorm(vmin=min(vmin, -1e-9), vcenter=0.0, vmax=max(vmax, 1e-9))


def _quiet_font_logging() -> None:
    """Silence fontTools' per-glyph INFO chatter when writing vector output.

    Writing a PDF emits several hundred INFO lines per figure describing font
    subsetting. Raised to WARNING rather than disabled, so genuine font errors
    (a missing family, a corrupt file) still surface.
    """
    for name in ("fontTools", "fontTools.subset", "fontTools.ttLib",
                 "matplotlib.font_manager"):
        logging.getLogger(name).setLevel(logging.WARNING)


def check_font_available() -> str | None:
    """Warn if neither Arial nor Helvetica is installed.

    DMM: "ARIAL or HELVETICA must be the font choice used throughout figure
    preparation." matplotlib falls through its sans-serif list silently, so on a
    machine without either font every figure renders in DejaVu Sans and looks
    fine locally while being non-compliant. Better to say so once, loudly.

    Returns the resolved family name, or None if it could not be determined.
    """
    from matplotlib.font_manager import FontProperties, findfont

    # Resolve against the configured family list rather than the generic
    # "sans-serif" alias, which findfont does not accept as a family name.
    try:
        family = mpl.rcParams.get("font.sans-serif") or list(COMPLIANT_FONTS)
        resolved = Path(findfont(FontProperties(family=family))).stem
    except Exception as exc:  # noqa: BLE001 - any failure here is cosmetic
        log.debug("Could not resolve the sans-serif font: %s", exc)
        return None

    if not any(f.lower() in resolved.lower() for f in COMPLIANT_FONTS):
        log.warning(
            "Figures will render in %r. DMM requires Arial or Helvetica; "
            "install one of them, or restyle the figures before submission.",
            resolved,
        )
    return resolved


def apply_style() -> None:
    """Set global matplotlib parameters to DMM figure specifications."""
    _quiet_font_logging()
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": DPI_SUBMISSION,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01,   # DMM: minimise background
        "savefig.facecolor": "white",
        "savefig.transparent": False,

        # DMM 4.3.1: 8 pt Arial for all labelling other than panel letters.
        # Sizes are set at the printed size, so nothing is scaled down later.
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": LABEL_PT,
        "axes.titlesize": LABEL_PT,
        "axes.titleweight": "bold",     # "headings should be bold"
        "axes.labelsize": LABEL_PT,
        "xtick.labelsize": LABEL_PT,
        "ytick.labelsize": LABEL_PT,
        "legend.fontsize": LABEL_PT,
        "legend.title_fontsize": LABEL_PT,
        "figure.titlesize": LABEL_PT,
        "figure.titleweight": "bold",

        "axes.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,

        "lines.linewidth": 0.75,
        "patch.linewidth": 0.4,
        "legend.frameon": False,

        # DMM: save text/font information as "text", not "curves" or "outlines",
        # so the journal can edit and reflow labels during production.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    check_font_available()
    log.debug("DMM figure style applied")


def fit_to_width(fig, width: float = DOUBLE_COL, tol: float = 0.002,
                 max_iter: int = 8) -> float:
    """Shrink a figure so its *rendered* width fits a journal column width.

    Needed because scanpy sizes multi-panel grids from the per-panel figsize and
    then adds legend and colourbar space outside the axes, so the final width is
    not knowable in advance. Left alone, a four-panel UMAP grid lands near
    385 mm against a 180 mm limit.

    Shrink-only: a figure already inside the target is returned untouched.
    Enlarging would spread a content-sized layout across empty page (DMM asks
    for the smallest size that conveys the information), and because fonts are
    fixed in points it also makes the type relatively smaller than designed.

    Iterative because font sizes are fixed in points: shrinking the canvas makes
    text occupy proportionally more space, which moves the content bbox, so a
    single scaling overshoots. Converges in two or three passes.

    Returns the achieved width in inches.
    """
    def rendered_width() -> float:
        fig.canvas.draw()
        try:
            return fig.get_tightbbox(fig.canvas.get_renderer()).width
        except Exception:  # noqa: BLE001 - renderer unavailable on some backends
            return fig.get_size_inches()[0]

    for _ in range(max_iter):
        current = rendered_width()
        if current <= 0:
            return fig.get_size_inches()[0]
        if current <= width * (1 + tol):
            break
        scale = width / current
        w, h = fig.get_size_inches()
        fig.set_size_inches(w * scale, h * scale)

    # The proportional iteration above does not always converge: for very wide,
    # short figures, shrinking the canvas makes the fixed-point type relatively
    # larger, which pushes the content bbox back out. Finish with a damped trim
    # that always undershoots, so the loop terminates from above rather than
    # oscillating around the target.
    final = rendered_width()
    for _ in range(4):
        if final <= width:
            break
        w, h = fig.get_size_inches()
        shrink = 0.99 * width / final
        fig.set_size_inches(w * shrink, h * shrink)
        final = rendered_width()

    if final > width:
        log.warning(
            "Figure could not be fitted below %.1f mm (stalled at %.1f mm); "
            "it likely has more panels or longer labels than fit the width. "
            "Split it or shorten the labels before submission.",
            width / MM, final / MM,
        )
    return final


def save_figure(
    fig_dir: Path,
    name: str,
    formats: tuple[str, ...] = DEFAULT_FORMATS,
    dpi: int = DPI_SUBMISSION,
    width: float | None = DOUBLE_COL,  # pass SINGLE_COL for half-page panels
    close: bool = True,
) -> list[Path]:
    """Write the current figure in each requested format at a journal width.

    PDF and SVG are the submission artefacts: DMM accepts EPS, PDF and SVG for
    line art and warns that JPEG or TIFF "may delay production". PNG is written
    alongside for review, preprints and README display only, and should not be
    what is uploaded at revision.
    """
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.gcf()
    if width is not None:
        assert_column_width(width, name)
        # savefig adds pad_inches on each side after the tight bbox is computed,
        # so fit to the target minus that padding or every figure lands ~1 mm over.
        pad = 2 * mpl.rcParams.get("savefig.pad_inches", 0.0)
        fit_to_width(fig, width - pad, tol=0.002, max_iter=8)
    else:
        # Content-sized figure (see dot_figsize): saved as designed, but the
        # page limit still applies, so measure rather than trust the estimate.
        fig.canvas.draw()
        try:
            rendered = fig.get_tightbbox(fig.canvas.get_renderer()).width
        except Exception:  # noqa: BLE001 - renderer unavailable on some backends
            rendered = fig.get_size_inches()[0]
        if rendered > MAX_WIDTH:
            log.warning(
                "%s renders at %.0f mm, over the %.0f mm maximum. It has more "
                "columns than fit the page at 8 pt; split it before submission.",
                name, rendered / MM, MAX_WIDTH / MM,
            )

    height_mm = fig.get_size_inches()[1] / MM
    if height_mm > MAX_HEIGHT / MM:
        log.warning("%s is %.0f mm tall, over the %.0f mm page maximum",
                    name, height_mm, MAX_HEIGHT / MM)
    elif height_mm > WORKING_HEIGHT / MM:
        log.info("%s is %.0f mm tall; the legend is set beneath the figure on "
                 "the page, so little room is left for it", name, height_mm)

    background = background_fraction(fig)
    if background > 0.55:
        log.info("%s is ~%.0f%% background. DMM asks that figures maximise data "
                 "and minimise white space; check for gaps between panels.",
                 name, 100 * background)

    written = []
    for fmt in formats:
        path = fig_dir / f"{name}.{fmt}"
        # metadata={} suppresses the software/creation-date tags matplotlib
        # would otherwise embed. DMM asks that no colour profile be assigned;
        # matplotlib writes plain RGB and attaches no ICC profile, so nothing
        # further is needed, but keeping the metadata empty also makes output
        # byte-reproducible across runs.
        kwargs = {"metadata": {}} if fmt in ("pdf", "svg", "png") else {}
        fig.savefig(path, dpi=dpi, format=fmt, **kwargs)
        written.append(path)
    if close:
        plt.close("all")
    return written


def scanpy_panel_size(n_cols: int = 2, total_width: float = DOUBLE_COL) -> None:
    """Size scanpy's auto-generated multi-panel grids to a journal column width.

    scanpy builds its own figure for sc.pl.umap and friends, sizing each panel
    from rcParams["figure.figsize"] rather than from any figure created here.
    Left at the matplotlib default this produces figures several hundred mm
    wide; scaled down for submission, the 7 pt labels set above become
    illegible. Setting the per-panel size so that `n_cols` panels fill the
    target width means the type is the size it claims to be.
    """
    panel = total_width / n_cols
    mpl.rcParams["figure.figsize"] = (panel, panel)


def assert_column_width(width: float, name: str = "figure") -> None:
    """Warn when a figure width is neither a single nor a double column.

    DMM's image guidelines state a preference for single- or double-column
    layouts and show a 1.5-column figure as a negative example, criticised for
    suboptimal data layout and unminimised background. Intermediate widths also
    force the typesetter to either scale the figure or leave dead page space.
    """
    if not any(abs(width - w) < 2 * MM for w in COLUMN_WIDTHS):
        log.warning(
            "%s is %.0f mm wide, which is neither a single (%.0f mm) nor a "
            "double (%.0f mm) column. DMM prefers one or the other.",
            name, width / MM, SINGLE_COL / MM, DOUBLE_COL / MM,
        )


def background_fraction(fig) -> float:
    """Estimate the proportion of the figure not covered by axes.

    A crude proxy for DMM's instruction to maximise information and minimise
    background within the figure area. Their negative examples are all
    criticised for white gaps between rows and ineffective use of space, so
    having a number makes that checkable rather than a matter of taste.

    Legends, colourbars and suptitles are drawn outside axes and count as
    background here, so a figure with a large external legend will score high
    even when well laid out. Treat it as a prompt to look, not a verdict.
    """
    fig.canvas.draw()
    fig_w, fig_h = fig.get_size_inches()
    total = fig_w * fig_h
    if total <= 0:
        return 0.0

    covered = 0.0
    for ax in fig.axes:
        bbox = ax.get_position()
        covered += (bbox.width * fig_w) * (bbox.height * fig_h)
    return max(0.0, 1.0 - covered / total)


def panel_spacing_for(adata, keys, base: float = 0.30, cap: float = 0.85) -> float:
    """Choose inter-panel spacing wide enough for the widest legend.

    scanpy draws each panel's legend in the gap to that panel's right, so the
    space a figure needs is set by its most category-rich and longest-labelled
    colour key, not by a number chosen once. A fixed value either clips a
    16-level cluster legend or wastes half the page on a two-level one.

    Scales with both the number of legend columns scanpy will use (it wraps
    after roughly eight entries) and the longest label, since a legend of
    "cmz_progenitor (low confidence)" is far wider than one of "day5".
    """
    widest = 0.0
    for key in keys:
        if key not in getattr(adata, "obs", {}):
            continue
        values = adata.obs[key]
        if not (hasattr(values, "cat") or values.dtype == object):
            continue
        levels = [str(v) for v in values.unique()]
        columns = max(1, -(-len(levels) // 8))          # scanpy wraps at ~8
        longest = max((len(v) for v in levels), default=0)
        # Per column: a fixed allowance for the marker and its gutter, plus a
        # term for the label text. Calibrated so that two columns of short
        # numeric cluster labels already exceed the base, since sixteen entries
        # need real width even when each is one character.
        widest = max(widest, columns * (0.16 + 0.011 * min(longest, 32)))

    return float(min(max(base, widest), cap))


def panel_label(ax, letter: str, x: float = -0.12, y: float = 1.06) -> None:
    """Add a panel letter in DMM house style: 12 pt bold uppercase."""
    ax.text(x, y, letter.upper(), transform=ax.transAxes,
            fontsize=PANEL_LETTER_PT, fontweight="bold", va="bottom", ha="right")


def figure(width: float = DOUBLE_COL, height: float | None = None, **kwargs):
    """Create a figure at a journal column width.

    Sizing at the final printed width means the font sizes set in `apply_style`
    are the sizes that appear in print. Drawing large and scaling down is the
    usual cause of unreadable 4 pt axis labels.
    """
    height = height or width * 0.6
    if width > MAX_WIDTH:
        log.warning("Figure width %.1f mm exceeds the DMM maximum of %.0f mm",
                    width / MM, MAX_WIDTH / MM)
    if height > MAX_HEIGHT:
        log.warning("Figure height %.1f mm exceeds the DMM maximum of %.0f mm",
                    height / MM, MAX_HEIGHT / MM)
    return plt.subplots(figsize=(width, height), **kwargs)
