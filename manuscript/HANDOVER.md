# What is real in this draft, and what is not

> **SUPERSEDED.** This is the pre-integration draft, written before the
> bulk RNA-seq, EM and IHC data were available. Use
> `results_discussion_integrated.md` and `PROVENANCE.md` instead.

Read this before using `scrnaseq_section.md`.

## Written from data

Everything in the Methods, the Results, and the first and last paragraphs of the
Discussion derives from the reanalysis in this repository. Numbers trace to
`results_cellbender_presupport/` and to `article/figS3_acceptance_criteria.csv`.

## Written as placeholders

The EM, IHC and bulk RNA-seq paragraphs are `[BRACKETED CAPITALS]`. **None of
those datasets were available to this analysis.** `bulk-rnaseq/` in the parent
directory is a workflow skeleton with no results, and no EM or IHC data is
present. The placeholders state what each paragraph must establish and how it
connects to the single-cell result; they are not summaries of findings, and
nothing should be inferred from their wording.

## Three things to resolve before submission

### 1. The acceptance-criteria table contains a false statement

`figS3_acceptance_criteria.csv` and `figS3_notes.txt` both say ambient-RNA
correction was "not applied in this run". **It was.** These outputs came from the
CellBender-corrected 22-library run. The text is hardcoded at `article.py:353`
rather than derived from the input, so it reports the same thing regardless.

Criterion 4 should read *pass*, and the S3 notes paragraph on ambient correction
should be deleted. Do not paste that table into the supplement as it stands — it
tells a reviewer the correction was not done.

### 2. Two criteria are marked FAIL and both are arguably mis-specified

**Criterion 1** ("each target detected in more RPE cells than *rho*") fails for
*slc39a8* (5.4% vs 8.3%) and *slc30a10* (2.2% vs 8.3%). *rho* is a poor sentinel
for this test: rods express it enormously, so its ambient floor is high in every
cluster and the bar is set arbitrarily high. Cluster-restricted enrichment
(criterion 3) is the better test and passes decisively — *slc39a8* is 75-fold
enriched in RPE, *slc30a10* is undetected outside it. Ambient contamination is
approximately uniform across cell types; a 75-fold cluster restriction is not
ambient. Consider replacing criterion 1 with an enrichment ratio against the
pooled non-RPE background.

**Criterion 6** ("reproduces at both stages") fails on *slc39a14*, 15.4% at 5 dpf
against 7.4% at 8 dpf. The 8 dpf libraries are shallower (median genes per cell
598–720 vs 686–863) and detection fraction scales with depth, so this may be
technical. A depth-matched comparison — downsampling the 5 dpf libraries to the
8 dpf median — would settle it. Until then the draft does not claim a
developmental change, and should not.

### 3. The figures predate the mixed-cluster check

Fig. 2D includes a `lens_fibre (low confidence)` row, and that cluster carries the
highest ambient-sentinel signal in the panel. The homogeneity check added on
2026-08-20 flags clusters whose cells mostly disagree with their own label; it had
not been run when these panels were generated. Re-running will not move the RPE
result — RPE clusters score 92–98% support — but it will relabel some rows in
Fig. 2D and Fig. S3A. Regenerate before submission so the panels and the code
agree.

## Figures as they stand

Publication-ready, from the ambient-corrected run, in
`results_cellbender_presupport/article/figures/` as PDF, SVG and PNG:

| Panel | File |
|---|---|
| 2C | `fig2C_rpe_by_stage`, `fig2C_rpe_markers` |
| 2D | `fig2D_transporter_controls` |
| S3A | `figS3_marker_dotplot_1` … `_4` |
| S3B | `figS3_cells_per_cluster` |
| S3C | `figS3_detection_fractions` |
| S3D | `figS3_acceptance_criteria.csv` — **fix criterion 4 first** |

Use the PDF or SVG for submission; the PNGs are for review copies.

## References

`references.ris` imports directly into Mendeley. Two entries are marked
`VERIFY page range` (Raj et al. 2020; Sur et al. 2023) — the DOIs are right and
Mendeley will correct the pagination on import, but check them.
