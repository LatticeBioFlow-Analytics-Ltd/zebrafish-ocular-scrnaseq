# screye — two-pass scRNA-seq analysis of the zebrafish head and eye

Two-pass QC, clustering and cell-type annotation of paired 10x libraries from
whole-head zebrafish dissociations at 5 dpf and 8 dpf. Pass 1 resolves tissue
compartments across the whole dissociation; pass 2 subsets the ocular and neural
compartments and re-clusters each at its own resolution.

Written as a small installable package with a config-driven CLI rather than a
notebook: the analysis is re-run whenever a sample or a threshold changes, and
the parameters need to live in one file rather than scattered through cells. A
walkthrough notebook in `notebooks/` calls the same functions for interactive
inspection.

**Data:** a reanalysis of Raj et al. 2020 (*Neuron* 108:1058; GEO GSE158142),
whole-head dissociations at two larval stages — see [`DATA.md`](DATA.md) for
provenance. Because the dissociation contains brain as well as eye, a retinal
marker panel alone cannot annotate it; that constraint is what the two-pass
design exists to handle.

## Run it

The environment is managed with [pixi](https://pixi.sh). Install it once:

```bash
curl -fsSL https://pixi.sh/install.sh | bash        # macOS / Linux
```

Then, from the repository root:

```bash
pixi run run -- --dry-run    # validate config and inputs, print the plan, run nothing
pixi run run                 # full pipeline on the real data (~35 min, see below)
pixi run -e dev test         # unit tests, no data required (~4 s)
pixi run -e dev integration  # full pipeline + notebook on synthetic data (~2 min)
pixi run -e dev check        # lint + unit tests
pixi task list               # everything available
```

There is no separate install or activate step: `pixi run` creates the
environment from `pixi.lock` on first use and reuses it afterwards. The lockfile
pins every package to an exact version and hash for linux-64, osx-64, osx-arm64
and win-64, so the environment is reproduced rather than re-solved — which is
the substantive difference from `conda env create -f environment.yml`, where two
people running the same command months apart get different builds.

The package installs itself into the environment as an editable dependency, so
`import screye` and the `screye` console script work without any `PYTHONPATH`
manipulation, including inside the notebook (`pixi run -e dev lab`).

`.pixi/config.toml` sets `detached-environments = true`, storing the
environments outside the repository. This repository lives in an
iCloud-synced folder, and iCloud both scatters `* 2.py`-style conflict copies
through `.pixi` and sets the macOS hidden flag on files mid-sync — Python
3.12.13+ deliberately skips hidden `.pth` files, which silently breaks the
editable install with `No module named 'screye'`. Detaching the environments
keeps them out of iCloud's reach (and out of your iCloud quota).

Everything else is driven by `config/config.yaml`. No threshold, path or seed
appears anywhere in the source.

<!-- 35 min measured 2026-08-12 on an Apple-silicon MacBook (16 GB), 53k cells
     in, t-SNE on ~50k cells is the dominant cost. Re-measure after any change
     to the embedding stage - a reviewer told "~3 min" who then waits 40 stops
     trusting the rest. -->

## How it fits together

Two passes. Pass 1 clusters the whole dissociation and assigns each cluster to a
tissue compartment; pass 2 subsets each compartment, re-selects features within
it, and re-clusters.

```
config/config.yaml ─────────┐
config/markers_compartment ─┤
config/markers_ocular ──────┤
                            v
data/*.h5 ──> qc.py ──> cluster.py ──> compartments.py ──> cluster.subcluster()
              QC,        whole-tissue    assign each         re-select HVGs,
              MAD        HVG, PCA,       cluster to a        re-embed, re-cluster,
              filters    UMAP, Leiden    compartment         fine annotation
                            │                 │                     │
                            v                 v                     v
                   results/figures/   pass1_compartment_    results/ocular/
                   (whole tissue)     evidence.csv          results/neural/
```

### Why two passes

Feature selection is a global operation. On a whole-head dissociation the top
2,000 highly variable genes are dominated by the axes separating brain from
muscle from blood; the genes distinguishing a bipolar cell from an amacrine cell
are not among them. Applying a retinal marker panel directly to this data does
not just fail to label non-ocular clusters, it forces them into whichever
retinal identity scores least badly.

Subsetting a compartment and re-running feature selection within it recovers the
structure that whole-tissue selection washes out. The same focused-clustering
approach is used in Sur et al. 2023 (Dev Cell 58:3028) and in Raj et al. 2020
(Neuron 108:1058), the source of these data.

Both compartments are retained. The `neural` subset is not a by-product: it is
the comparison that makes the ocular result interpretable, since a gene that
looks ocular-enriched against whole-tissue background may simply be
neuron-enriched.

### Marker panels

| Panel | Pass | Purpose |
|---|---|---|
| `config/markers_compartment.yaml` | 1 | 15 broad tissue compartments — ocular, neural, immune, muscle, blood, cartilage, epidermis |
| `config/markers_ocular.yaml` | 2 | 17 ocular cell types — photoreceptor classes, retinal interneurons, RGC, Müller glia, RPE, CMZ progenitors, lens, cornea |

Primary sources for the retinal sets: `rlbp1a`/`apoeb` for Müller glia
(Bernardos et al. 2007, J Neurosci); `vsx1` for bipolar cells (Vitorino et al.
2009, Neural Dev); `gad1b`/`gad2` for amacrine cells (Sandell et al. 1994);
`pde6` family for photoreceptors (Abalo et al. 2020); `rbpms2b`/`isl2b` for
retinal ganglion cells (Kölsch et al. 2021, Neuron 109:645); `cldn19` for
endothelium (Kolosov et al. 2013). Broader context from Hoang et al. 2020
(Science 370:eabb8598) and Daniocell (Sur et al. 2023).

Both panels are literature-assembled starting points, **not validated
references**. Every symbol should be checked against ZFIN and against Daniocell
at the matched stage before an annotation is relied on. Kölsch et al. resolved
>30 RGC types in zebrafish; this pipeline deliberately does not attempt that —
at 5–8 dpf with these cell numbers, class-level resolution is the honest ceiling.

## Article figure panels (Fig 2C, 2D, S3)

The manuscript panels are generated from the pass-2 ocular subset, either as
part of a full run (via `article_panel:` in the config) or standalone in
seconds from the saved subset:

```bash
pixi run -e dev python -m screye.article --xlsx path/to/transporter_dotplot_control_panel.xlsx
```

Outputs land in `results/article/`: the RPE UMAP per stage with the cluster
size stated (2C), the 12-gene control-panel dotplot with dot size as detection
fraction on an absolute scale and colour as mean expression among expressing
cells only (2D), and the S3 set — cluster markers, cells per cluster per
stage, RPE vs non-RPE detection fractions — plus `figS3_acceptance_criteria.csv`,
an auto-evaluated version of the control workbook's criteria, and (with
`--xlsx`) a filled copy of that workbook. The gene panel and its tiers live in
`config/article_panel.yaml`, transcribed from the workbook; two criteria
(ambient-RNA correction, second-atlas replication) need external inputs and
are recorded as pending in `figS3_notes.txt`.

## Data layout

```
data/
  day5_filtered_feature_bc_matrix.h5
  day8_filtered_feature_bc_matrix.h5
```

The mapping from filename to sample name and timepoint is explicit in
`config/config.yaml` rather than parsed from the filename, so renaming a file
cannot silently relabel a condition.

### Consuming the `scRNA_analysis` Snakemake workflow

The upstream half of this project — SRA download through Cell Ranger and
CellBender, in the sibling `scRNA_analysis` repository — publishes three input
tiers under `scRNA_analysis/results/h5/`: per-library ambient-corrected
matrices (`cellbender/`), per-library raw matrices (`raw/`), and per-timepoint
aggregates. `config/config_snakemake.yaml` points this pipeline at the
CellBender tier, which is that workflow's stated analysis input:

```bash
pixi run run -- -c config/config_snakemake.yaml --dry-run
```

Instead of listing 22 libraries by hand, the config names the Snakemake
project's own `config/libraries.tsv` as a `sample_sheet:` plus an `h5_pattern:`
template — the sheet's `library` and `timepoint` columns stay the single source
of truth, and nothing is parsed from filenames. Every figure scales with the
library count: violins widen, the composition plot shows one point per library
(real spread, now that there is one), and stacked bars shade libraries by
timepoint. With replication per timepoint, `sample` is no longer confounded
with stage, so `cluster.batch_correct` becomes a legitimate option there —
still off by default.

## Figures

Written to the Disease Models & Mechanisms figure specification
([manuscript preparation](https://journals.biologists.com/dmm/pages/manuscript-prep),
section 4, last updated 22 July 2026):

| Requirement | Implementation |
|---|---|
| Maximum 180 mm × 210 mm including lettering | Every figure fitted to 179.5 mm; both ceilings asserted in tests |
| Single- or double-column layouts preferred | Widths constrained to 85 mm or 180 mm; anything intermediate warns, since DMM shows a 1.5-column figure as a negative example |
| Maximise data, minimise background | Panel spacing tightened throughout; each save reports the background fraction and flags figures above 55% |
| Legible at printed size | Grid figures (dotplots, heatmaps) are sized from their row and column counts at a fixed ~3.3 mm pitch and saved as designed, never rescaled: fonts and dots are fixed in points, so shrinking an oversized canvas is what makes labels and dots collide. Panels wider than the 180 mm budget are split into numbered figures with all clusters kept in every part |
| Readable cluster keys | Cluster panels label clusters on the data; long categorical keys (cell types, compartments) get a panel with reserved legend space instead of a legend drawn over the neighbouring panel |
| 12 pt bold uppercase panel letters | `figures.panel_label()` |
| 8 pt Arial labelling, bold headings | Set in `figures.apply_style()` |
| Arial or Helvetica throughout | Checked at style time; warns and names the resolved font if neither is installed |
| Line art as EPS, PDF or SVG | PDF and SVG written for every figure. PNG is for review and preprints only — DMM warns that JPEG/TIFF for line art may delay production |
| Text saved as text, not outlines | `svg.fonttype: none`, `pdf.fonttype: 42`; asserted in tests |
| RGB, no spot colours, no colour profile | matplotlib writes plain RGB; no ICC profile attached |
| Show the true data spread | Composition plots show per-library points with cell counts rather than summary bars; violin box-and-whisker definitions written to a `*_legend.txt` beside each figure |

**Colour convention.** Diverging quantities use **red for positive/increased and
blue for negative/decreased**, centred on zero via `TwoSlopeNorm` so white always
means "no change". Without the centring, a −1 to +4 range places white at +1.5
and renders every genuinely unchanged feature pale blue — a visual claim of
decrease the data do not support.

This is chosen to satisfy DMM's colour-vision guidance, which asks authors to
avoid red/green for two-channel displays and encourages colour-blind-safe
choices particularly in heatmaps. Red-to-blue is the recommended alternative
axis rather than the prohibited one; the palette is ColorBrewer RdBu reversed.
Sequential quantities use viridis; categorical variables use Okabe-Ito.

**Sizing.** Figures are designed at their printed size rather than drawn large
and rescaled, because fonts and scatter dots are fixed in points: shrinking a
385 mm scanpy grid to 180 mm leaves the 8 pt type and the dots at full size
while the space between them halves, which is how dots come to overlap and
legends to collide. Embedding figures are laid out at column width directly;
dotplots and heatmaps get a canvas computed from their row and column counts at
a fixed ~3.3 mm pitch and are saved as designed. A damped iterative fit
(`figures.fit_to_width`) remains only as a trim for figures already close to
target, since a naive proportional loop oscillates rather than converging.

**Legend space.** DMM's page schematics place the figure legend in its own box
beneath the figure, so a figure using the full 210 mm leaves no room for it.
`WORKING_HEIGHT = 175 mm` is the height to design to; 210 mm remains the hard
ceiling and is enforced.

**Verify before submission.** The specification above is transcribed from the
[manuscript preparation guidelines](https://journals.biologists.com/dmm/pages/manuscript-prep)
and the [image guidelines PDF](https://journals.biologists.com/DocumentLibrary/DMM/DMM%20image%20guidelines.pdf)
as of August 2026. Journal requirements change; re-check both before submitting.

**Image manipulation.** DMM screens all accepted manuscripts for signs of image
manipulation and can reject on non-compliance. Nothing in this pipeline
manipulates micrographs, but if figures from here are assembled alongside
photographic panels: adjustments must be applied to the whole image, non-linear
adjustments such as gamma changes must be disclosed in both the figure legend
and Materials and Methods, and any cropping or consolidation must be visibly
marked and stated in the legend.

## Design decisions

**Batch correction is off by default.** The two samples *are* the two
developmental stages, so `sample` and `timepoint` are the same variable.
Correcting on `sample` would remove exactly the signal the experiment measures.
The option exists (`cluster.batch_correct`) for a genuine technical batch, and
logs a warning when used. `cluster_composition_by_sample.png` makes the
alternative reading checkable: a cluster drawn overwhelmingly from one library is
either stage-specific biology or an uncorrected artefact, and those need telling
apart before any label is trusted.

**Mitochondrial genes use the lower-case `mt-` prefix.** Zebrafish symbols follow
ZFIN convention, so the human `MT-` prefix matches nothing and returns a
mitochondrial fraction of exactly zero for every cell — a filter that appears to
run and silently does not. `annotate_gene_classes` warns when no genes match, and
a unit test pins the behaviour in both directions.

**QC thresholds are in median absolute deviations, not fixed counts.** Fixed
cut-offs copied between datasets are the usual cause of quietly discarding a real
cell type with low RNA content, and the two libraries here differ in depth. A
hard 20% mitochondrial ceiling sits on top as a backstop.

**Cells are flagged before they are filtered.** `qc_fail` is written to `.obs`
and the QC figures are drawn pre-removal, so it is visible when a whole
population is about to be dropped. The manifest records per-sample attrition.

**Doublets are scored per sample, on QC-passing cells only.** Pooling libraries
first would let barcodes from one library act as simulated doublet partners for
the other. Restricting to QC-passing cells keeps near-empty droplets out of the
simulated pairs — and is required on ambient-corrected input, where a barcode
can be stripped below Scrublet's own minimum-count filter.

**The doublet threshold is checked, not trusted.** Scrublet picks one with
Otsu's method on the simulated-doublet score histogram, which presumes that
histogram is bimodal. On shallow data from a continuous tissue it is not, and
the method can return a threshold above every observed cell — calling nothing,
indistinguishable in the output from a dataset containing no doublets. That is
refused rather than reported: set `qc.doublet_threshold` to a value read off
the histogram, or set `qc.run_doublet_detection: false` and say so. Whichever
happens is recorded in the run manifest, so a run that skipped the step cannot
be mistaken for one that found nothing.

**Leiden runs at three resolutions,** all retained in `.obs`, so the choice can
be justified against the data rather than asserted.

**Annotation is auditable, not assertive.** Each cluster is scored against every
marker set and labelled by its top-scoring set, but the full score matrix is
written out, and any cluster whose top two scores differ by less than 0.05 is
labelled `(low confidence)`. A cluster scoring nearly equally for two types is
one the panel cannot resolve; saying so is more useful than a confident wrong
label.

**The marker panel includes non-retinal populations** — lens, cornea, periocular
mesenchyme, erythrocytes, skeletal muscle. A whole-eye dissociation is not a
retinal dissociation, and omitting those sets does not remove the cells, it
forces them into whichever retinal label scores least badly.

**Both UMAP and t-SNE are produced** because they disagree usefully: t-SNE tends
to fragment continuous trajectories that UMAP joins. Neither preserves global
inter-cluster distance, so proximity in either is not evidence of relatedness,
and neither supports any claim here beyond visual inspection.

## Limitations

- `config/markers.yaml` is a literature-assembled starting point, not a validated
  reference. Every symbol should be checked against ZFIN and against the
  Daniocell atlas (Sur et al., 2023) at the matched stage.
- Ambient RNA is not corrected. For whole-eye dissociations, CellBender or SoupX
  would be the next addition and would likely reduce apparent photoreceptor
  transcript expression in non-photoreceptor clusters.
- Marker-based labelling is coarse. Reference mapping onto Daniocell with scVI or
  `sc.tl.ingest` would resolve closely related interneuron subtypes better, at
  the cost of inheriting that atlas's own annotation errors.
- One library per timepoint means no biological replication: compositional
  differences between stages are descriptive, not tested.
