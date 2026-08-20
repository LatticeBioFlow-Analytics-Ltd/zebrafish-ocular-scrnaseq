# Where every number in the draft comes from

## Single-cell (Fig. 2C–D, S3)

`~/Desktop/Lattice BioFlow/screye/results_cellbender/`

| Claim | File |
|---|---|
| 107,444 cells, 34 clusters | `run_manifest.json` |
| ocular subset 4,572 cells | `run_manifest.json` → `pass2_subsets.ocular` |
| RPE 246 / 163 cells | `ocular/celltype_counts.csv`, summed by stage |
| detection fractions, all controls | `article/figS3_acceptance_criteria.csv` |
| reads-in-cells 0.49–0.56 / 0.33–0.47 | upstream Cell Ranger `metrics_summary.csv` |

## Bulk RNA-seq (Fig. 2A–B, 3, 4)

`~/Desktop/Bioinformatics/EyeData/FINAL_RESULTS/`

| Claim | File |
|---|---|
| DEG counts, all four contrasts | `dge/DEG_counts_summary.csv` |
| n = 3 per group, fully crossed | `dge/QC_replicate_by_group.csv` |
| *atf3*, *atf5a*, *jun*, *creb5a*, *pde6h*, *grk1b*, *lrit1b* | `dge/DE_interaction_GxT.csv` |
| transporters not DE (P adj > 0.4) | `dge/DE_*.csv`, all four |
| GO NES and P adj | `gsea/interaction_GO_GSEA.csv` |
| KEGG NES and P adj | `kegg/interaction_KEGG_GSEA.csv` |
| cell-type log FC | `celltype/celltype_diff_GSVA.csv`, `contrast == "interaction"` |

## EM (Fig. 6)

`~/Desktop/zf-em-mito-annotation/data/results/`

| Claim | File |
|---|---|
| group means, area fraction 0.195 → 0.061, 0.116 → 0.200 | `group_summary.csv` |
| all Hedges' *g* and CIs | `effect_sizes.csv` |
| 25,000× the only comparable magnification | `README.md`, magnification × group table |
| MitoNet artefact, *g* = +2.9, ρ = +0.34 / +0.53 | `README.md` |

## IHC — **no measurements exist**

`~/Desktop/Bioinformatics/EyeData/Manuscript/ihc_analysis_v2/`

240 images / 239 fish, one image per fish, blinded (`BLINDING_KEY.csv`),
`wt` 114 / `mut` 126 across three rounds. Marker panels from `manifest_clean.csv`:

```
GNAT2_4C4 55   RIBEYE_GS 50   CASP3_ZRF1 39   RPE65_HUCD 35
PKC_HUCD  30   RPE65_ZRF1 16  CASP3_ZPR2 15
```

`03_measure.ijm` and `04_stats_superplots.R` exist but have produced no output.
The draft therefore reports the IHC as **in progress** and makes no claim from
it. Do not upgrade that sentence without the measurements.

---

# Things to fix before submission

## 1. RESOLVED — the acceptance-criteria table is now correct

Criterion 4 read "not applied in this run" when the run *was* CellBender-
corrected. The string was hardcoded in two places (`_acceptance_criteria` and
the `figS3_notes.txt` block); both now derive from `qc.ambient_corrected`,
declared in the config and settable on the standalone CLI via
`--ambient-corrected`. It is deliberately not inferred from the input path.
Criterion 4 now reads **pass**, and the notes paragraph reads APPLIED.

Panels and tables were regenerated on 2026-08-20 from
`results_cellbender/ocular/ocular.h5ad`.

## 2. RESOLVED — the two FAIL verdicts now carry their reasoning

The verdicts are unchanged: the control workbook defines the criteria and this
pipeline does not redefine them. But every row now carries a computed `note`, so
a bare FAIL in a supplement is interpretable.

**Criterion 1** records that *rho* is a conservative sentinel (rods express it at
extreme level, raising its ambient floor everywhere) and gives the enrichment
ratios that criterion 3 measures: *slc39a14* 3x, *slc39a8* **75x**, *slc30a10*
undetected outside the RPE.

**Criterion 6** records that it inherits the *rho* comparison, and reports the
per-stage depth measured from the data — median genes per cell 622 at 5 dpf
against 537 at 8 dpf, a 1.16x difference — so the stage difference is not read
as developmental.

Note the enrichment for *slc39a8* is **75x**, not the 54x quoted in the first
draft: 54 came from dividing the rounded display values (5.4 / 0.1), whereas the
unrounded figures are 5.379% and 0.072%.

## 3. MuSiC deconvolution fits are poor

`celltype/music/music_fit_quality.csv` reports per-sample R² of 0.18–0.20. The
draft reports the **GSVA** signature scores, which do not depend on that fit. If
MuSiC proportions appear in any figure, the R² should be stated or the panel
dropped.

## 4. EM has no animal-level identifier

`data/manifest/frames.csv` indexes frames by group directory only
(`2024-10-16 SLC39a14 6dpf MnCl2 exposed 22` etc.). Fish identity is not
recorded, so *n* is fields, not animals, and within-animal correlation cannot be
excluded. The draft says so explicitly. If fish can be recovered from acquisition
metadata or lab notes, a mixed model with fish as a random effect would
strengthen Fig. 6 considerably.

## 5. RESOLVED — figures carry the mixed-cluster flags

The pipeline was re-run on 2026-08-20 with the homogeneity check active. Nothing
about the analysis moved: 100% of cells stayed in the same Leiden cluster,
107,444 cells and 34 clusters as before, identical compartment split, and every
relabelling was `X` to `X (mixed)`. The acceptance-criteria table and detection
fractions are byte-identical to the previous run.

The RPE cluster is unflagged (92-98% support) and unchanged at 409 cells, so
Fig. 2C and 2D are unaffected. `results_cellbender_presupport/` was kept only
long enough to prove that and has been removed.

One finding came out of the re-run and is worth carrying into the Methods: the
whole-head cell-type labels are scored against the *ocular* panel, and 16 of 34
of those clusters fall below 50% support, against 4 of 34 for the compartment
panel on identical clusters. Fine cell types are only meaningful within a
compartment; the trustworthy whole-head annotation is the compartment label.
