"""Typed configuration objects loaded from YAML.

All analysis parameters live in `config/config.yaml` so that no threshold,
path or random seed is ever hard-coded in the analysis modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class SampleSpec:
    """One CellRanger output to be loaded.

    Attributes
    ----------
    name:
        Short sample identifier used in plots and `adata.obs["sample"]`.
    h5:
        Path to the CellRanger `filtered_feature_bc_matrix.h5`.
    timepoint:
        Biological stage label, e.g. ``"5dpf"``. Kept separate from `name`
        because timepoint is a biological covariate, not a technical batch.
    """

    name: str
    h5: Path
    timepoint: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "h5", Path(self.h5))


@dataclass(frozen=True)
class SampleSheet:
    """Samples declared by reference to an upstream sample sheet.

    This is the handover from the `scRNA_analysis` Snakemake workflow, whose
    `config/libraries.tsv` is the source of truth for which 10x libraries exist
    and which timepoint each belongs to. Re-declaring 22 libraries by hand in
    `samples:` invites exactly the copy-paste drift the explicit mapping was
    meant to prevent, so the sheet is read instead and the h5 path derived from
    a template. Nothing is parsed out of filenames: the timepoint comes from
    the sheet's own `timepoint` column.

    Attributes
    ----------
    libraries:
        TSV with at least `library` and `timepoint` columns (the Snakemake
        project's `config/libraries.tsv` as-is).
    h5_pattern:
        Path template with a ``{library}`` placeholder, e.g.
        ``../scRNA_analysis/results/h5/cellbender/{library}_cellbender_filtered.h5``.
    """

    libraries: Path
    h5_pattern: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "libraries", Path(self.libraries))
        if "{library}" not in self.h5_pattern:
            raise ValueError(
                "sample_sheet.h5_pattern must contain a {library} placeholder, "
                f"got {self.h5_pattern!r}"
            )

    def expand(self, root: Path) -> list[SampleSpec]:
        """Read the sheet and return one SampleSpec per library, in sheet order."""
        import csv

        sheet = self.libraries if self.libraries.is_absolute() else root / self.libraries
        if not sheet.exists():
            raise FileNotFoundError(f"sample sheet not found: {sheet}")

        with sheet.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            missing = {"library", "timepoint"} - set(fields)
            if missing:
                raise ValueError(
                    f"{sheet}: missing column(s) {sorted(missing)}; "
                    f"found {fields}"
                )
            rows = [r for r in reader if r.get("library")]

        if not rows:
            raise ValueError(f"{sheet}: no library rows")

        specs = []
        for row in rows:
            h5 = Path(self.h5_pattern.format(library=row["library"]))
            specs.append(SampleSpec(
                name=row["library"],
                h5=h5 if h5.is_absolute() else root / h5,
                timepoint=row["timepoint"],
            ))
        return specs


@dataclass(frozen=True)
class QCConfig:
    """Quality-control parameters.

    Thresholds are expressed in median absolute deviations (MADs) rather than
    fixed counts. Fixed cut-offs copied between datasets are the most common
    source of silently discarded cell types; MAD-based cut-offs adapt to each
    sample's own distribution.
    """

    mito_prefix: str = "mt-"  # zebrafish gene symbols are lower-case
    ribo_prefixes: tuple[str, ...] = ("rps", "rpl")
    n_mads: float = 5.0
    max_pct_mito: float = 20.0  # absolute ceiling applied on top of the MAD rule
    min_genes_per_cell: int = 200
    min_cells_per_gene: int = 3
    run_doublet_detection: bool = True
    expected_doublet_rate: float = 0.06
    # Score above which a cell is called a doublet; "auto" uses Scrublet's own
    # threshold. That automatic choice is Otsu's method applied to the
    # simulated-doublet score histogram, which presumes the histogram is
    # bimodal — a low "embedded" mode of doublets indistinguishable from
    # singlets, and a high "neotypic" mode that is detectable.
    #
    # On shallow data from a continuous tissue the histogram is unimodal with a
    # long tail, and Otsu has no trough to find. On GSE158142 it returned 0.825
    # for both timepoints: above the highest-scoring observed cell (0.729) and
    # above 100% of the simulated doublets it was derived from, so no cell could
    # ever be called and none was. A threshold excluding every one of the
    # method's own positive controls is degenerate, and `detect_doublets`
    # refuses it rather than reporting zero doublets as a finding.
    doublet_threshold: float | str = "auto"


@dataclass(frozen=True)
class ClusterConfig:
    """Dimensionality reduction and clustering parameters."""

    n_top_genes: int = 2000
    n_pcs: int = 30
    n_neighbors: int = 15
    leiden_resolutions: tuple[float, ...] = (0.4, 0.8, 1.2)
    primary_resolution: float = 0.8
    compute_tsne: bool = True
    tsne_perplexity: float = 30.0
    # Batch correction is OFF by default. See README: the two samples differ by
    # developmental stage, so "batch" and "biology" are confounded and blind
    # correction would remove the signal of interest.
    batch_correct: bool = False
    batch_key: str = "sample"


@dataclass(frozen=True)
class CompartmentConfig:
    """Two-pass compartment assignment and sub-clustering parameters."""

    # Minimum gap between the best and second-best compartment score for an
    # assignment to be treated as confident. Below this the cluster keeps the
    # top label but is flagged, and the flag propagates into the outputs.
    min_margin: float = 0.05
    # Minimum fraction of a cluster's cells that must individually score highest
    # for the label the cluster was given. Below this the cluster is a mixture
    # and the label is suffixed `(mixed)`.
    #
    # This is a different question from the margin, and the margin cannot answer
    # it: margins compare panels against each other on a cluster mean, so a
    # cluster containing two distinct populations can score one of them
    # decisively and never be flagged. In the uncorrected GSE158142 ocular
    # subset, a cluster of 853 cells holding 335 cones and 129 lens fibre cells
    # was labelled `lens_fibre` on a margin of 0.152 — three times the threshold
    # above — because crystallins are 35-40% of a lens cell's UMIs and lift the
    # mean on their own. Support for that cluster is ~0.28.
    #
    # 0.5 is a deliberately permissive floor: it asks only that a plurality of
    # cells agree with their own label. Raise it for a homogeneous tissue; a
    # continuum of related progenitor states will legitimately sit lower.
    min_support: float = 0.5
    # Compartments to subset and re-cluster in pass 2. Both are kept: the
    # neural subset is the comparison that makes the ocular result meaningful.
    subset_groups: tuple[str, ...] = ("ocular", "neural")
    # Below this many cells, sub-clustering describes sampling noise.
    min_cells_for_subset: int = 200
    # Pass-2 clustering is finer than pass 1: within a single compartment the
    # relevant distinctions are between related cell types, not between tissues.
    leiden_resolutions: tuple[float, ...] = (0.4, 0.8, 1.2)
    primary_resolution: float = 0.8
    n_top_genes: int = 2000
    n_pcs: int = 30
    n_neighbors: int = 15


@dataclass(frozen=True)
class Config:
    """Top-level configuration."""

    samples: list[SampleSpec]
    outdir: Path
    markers_file: Path            # pass-2 ocular panel (fine cell types)
    compartment_markers: Path     # pass-1 compartment panel (broad tissues)
    # Alternative to `samples`: expanded into SampleSpecs by `resolve()`, since
    # reading the sheet needs a known root to resolve relative paths against.
    sample_sheet: SampleSheet | None = None
    # Optional: gene panel for the article figure set (Fig 2C/2D/S3). When set,
    # the pipeline builds those panels from the ocular subset after pass 2.
    article_panel: Path | None = None
    seed: int = 0
    # Genes to plot individually across clusters and cell types, e.g. the
    # transporters whose expression pattern motivated the experiment.
    genes_of_interest: tuple[str, ...] = ()
    qc: QCConfig = field(default_factory=QCConfig)
    cluster: ClusterConfig = field(default_factory=ClusterConfig)
    compartment: CompartmentConfig = field(default_factory=CompartmentConfig)

    @classmethod
    def from_yaml(cls, path: str | Path) -> Config:
        """Load and validate a configuration file."""
        path = Path(path)
        with path.open() as handle:
            raw: dict[str, Any] = yaml.safe_load(handle)

        has_samples = bool(raw.get("samples"))
        has_sheet = bool(raw.get("sample_sheet"))
        if has_samples == has_sheet:
            raise ValueError(
                f"{path}: define exactly one of `samples` (explicit list) or "
                "`sample_sheet` (expanded from an upstream libraries TSV)"
            )

        samples = [SampleSpec(**s) for s in raw.get("samples") or []]
        sheet = SampleSheet(**raw["sample_sheet"]) if has_sheet else None
        qc = QCConfig(**raw.get("qc", {}))
        clustering = raw.get("cluster", {})
        for key in ("leiden_resolutions",):
            if key in clustering:
                clustering[key] = tuple(clustering[key])
        cluster = ClusterConfig(**clustering)

        compartment_raw = raw.get("compartment", {})
        for key in ("leiden_resolutions", "subset_groups"):
            if key in compartment_raw:
                compartment_raw[key] = tuple(compartment_raw[key])
        compartment = CompartmentConfig(**compartment_raw)

        return cls(
            samples=samples,
            sample_sheet=sheet,
            outdir=Path(raw["outdir"]),
            markers_file=Path(raw["markers_file"]),
            compartment_markers=Path(raw["compartment_markers"]),
            article_panel=(Path(raw["article_panel"])
                           if raw.get("article_panel") else None),
            seed=int(raw.get("seed", 0)),
            genes_of_interest=tuple(raw.get("genes_of_interest", ())),
            qc=qc,
            cluster=cluster,
            compartment=compartment,
        )

    def resolve(self, root: str | Path) -> Config:
        """Return a copy with every path made absolute against `root`.

        Paths in `config/config.yaml` are written relative to the repository
        root because that is where the CLI is run from. A notebook has a
        different working directory, and an absolute path in the YAML would not
        survive being shared. Resolving explicitly against a known root means
        one config file serves both without either carrying machine-specific
        paths.
        """
        root = Path(root).resolve()

        def _abs(path: Path) -> Path:
            path = Path(path)
            return path if path.is_absolute() else root / path

        # A sample sheet expands here rather than in `from_yaml` because the
        # sheet path is repository-relative and only the caller knows the root.
        samples = (
            self.sample_sheet.expand(root) if self.sample_sheet
            else [replace(s, h5=_abs(s.h5)) for s in self.samples]
        )

        return replace(
            self,
            samples=samples,
            outdir=_abs(self.outdir),
            markers_file=_abs(self.markers_file),
            compartment_markers=_abs(self.compartment_markers),
            article_panel=_abs(self.article_panel) if self.article_panel else None,
        )

    def validate_inputs(self) -> None:
        """Fail early and loudly if an input file is missing."""
        if self.sample_sheet and not self.samples:
            raise ValueError(
                "sample_sheet is declared but not expanded; call "
                "Config.resolve(root) before validate_inputs()/run()."
            )
        missing = [s.h5 for s in self.samples if not s.h5.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing CellRanger .h5 file(s):\n  "
                + "\n  ".join(str(m) for m in missing)
                + "\n\nRun `python tests/fixtures/make_synthetic_h5.py` to generate synthetic "
                "synthetic dataset, or edit config/config.yaml to point at your own."
            )
