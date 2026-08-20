# scRNA-seq reanalysis — draft text for DMM

> **SUPERSEDED.** This is the pre-integration draft, written before the
> bulk RNA-seq, EM and IHC data were available. Use
> `results_discussion_integrated.md` and `PROVENANCE.md` instead.

Draft of the Methods, Results and Discussion covering Fig. 2C–D and Fig. S3.
Everything here derives from the reanalysis in this repository. Passages in
`[SQUARE BRACKETS AND CAPITALS]` are placeholders for the EM, IHC and bulk
RNA-seq results, which were not available to the analysis and must not be
written from assumption — see `HANDOVER.md`.

---

## Methods

### Single-cell RNA-sequencing reanalysis

Publicly available single-cell RNA-sequencing data from dissociated zebrafish
heads at 5 and 8 days post-fertilisation (dpf) were obtained from the Sequence
Read Archive (BioProject PRJNA664124; GEO GSE158142) (Raj et al., 2020). The 136
sequencing runs correspond to 22 10x Genomics Chromium 3′ v2 libraries, each
sequenced four to eight times. Re-sequencing runs of a single library were
combined within one `cellranger count` invocation rather than aggregated
afterwards, because unique molecular identifiers can only be collapsed within a
single execution of the algorithm; aggregating per run inflates the recovered
cell number approximately 2.75-fold.

Reads were aligned and counted with Cell Ranger 9.0.1 (Zheng et al., 2017)
against GRCz11 (Ensembl release 116) filtered to protein-coding genes (25,432
genes, including all 13 mitochondrial genes). Ambient RNA was removed per library
with CellBender `remove-background` v0.3.2 (Fleming et al., 2023), using Cell
Ranger's own cell call as the expected cell number. Correction was applied
because the fraction of reads in cells differed systematically between stages
(0.49–0.56 at 5 dpf; 0.33–0.47 at 8 dpf, non-overlapping ranges), so ambient
contamination would otherwise be confounded with developmental stage.

Downstream analysis used Scanpy 1.12 (Wolf et al., 2018). Cells were filtered on
median-absolute-deviation outlier thresholds for total counts, genes detected and
mitochondrial fraction (5 MADs), with an absolute ceiling of 20% mitochondrial
reads and a minimum of 200 genes per cell (Luecken and Theis, 2019). Doublets
were scored per library with Scrublet (Wolock et al., 2019). Scrublet's automatic
threshold was unusable on these data — the simulated-doublet score histogram is
unimodal, so the threshold returned (0.825) exceeded every observed cell score
(maximum 0.729) and called nothing — and an explicit threshold of 0.25 was used
instead, selected as the point at which the observed call rate agrees with the
rate predicted from the assumed doublet rate and the measured sensitivity. This
removes the detectable minority of doublets (880 cells, 0.81%, at 12%
sensitivity); undetectable doublets remain and are reported as a limitation.

Of 208,308 barcodes, 107,444 cells were retained. Counts were library-size
normalised to 10,000 and log-transformed, 2,000 highly variable genes selected,
and Leiden clustering (Traag et al., 2019) applied at resolution 0.8 over 30
principal components, visualised with UMAP (McInnes et al., 2018). Clusters were
assigned to broad tissue compartments and then, within the ocular compartment,
re-analysed from feature selection onward and annotated against a curated marker
panel using `score_genes` (Satija et al., 2015). Batch correction was not
applied. Every cluster label carries two recorded diagnostics: the margin between
the best and second-best marker-set score, and the fraction of the cluster's
cells that individually score highest for the assigned label.

Detection fraction — the proportion of cells with at least one count for a gene —
is reported in preference to mean expression, because for low-abundance
transcripts a mean taken across all cells conflates abundance with dropout.

Analysis code, configuration and per-run provenance manifests are available at
[REPOSITORY URL / DOI].

---

## Results

### Manganese transporters are enriched in the retinal pigment epithelium

To ask whether the ocular manganese accumulation we observe [CROSS-REF FIG. 2A–B]
is accompanied by local manganese-handling capacity, we reanalysed a public
single-cell atlas of the zebrafish head at 5 and 8 dpf (Raj et al., 2020),
resolving 107,444 cells into 34 clusters and sub-clustering the ocular
compartment (4,572 cells; Fig. 2C, Fig. S3).

A single population expressing *rpe65a*, *tyrp1b* and *pmela* was identified as
retinal pigment epithelium (RPE) at both stages (246 cells at 5 dpf; 163 at
8 dpf; Fig. 2C). All three manganese transporters were detected in this
population and were enriched relative to every other ocular cell type: *slc39a14*
in 12.2% of RPE cells against 3.7% elsewhere, *slc39a8* in 5.4% against 0.1%, and
*slc30a10* in 2.2% against 0.0% (Fig. 2D). *slc39a8* and *slc30a10* were
essentially confined to the RPE. Detection of *slc39a14* was comparable to the
abundance-matched positive control *slc11a2* (12.5%), while out-of-tissue
controls remained near zero in every cluster (*cryaa* ≤3.3%; *hbae1.1* ≤0.6%),
indicating that the signal is not an artefact of residual ambient RNA (Fig. 2D,
Fig. S3).

Two features of the data limit what can be concluded. First, absolute detection
of *slc39a8* and *slc30a10* falls below that of the ambient sentinel *rho* within
the RPE (8.3%); because *rho* is expressed at very high levels in rods, its
ambient floor is high in every cluster, and cluster-restricted enrichment — 75-fold
for *slc39a8* and complete for *slc30a10* — is the more informative test, since
ambient contamination is approximately uniform across cell types. Second,
detection of *slc39a14* is lower at 8 dpf than at 5 dpf (7.4% versus 15.4%), but
the 8 dpf libraries are also shallower (median genes per cell 598–720 versus
686–863), and detection fraction scales with depth; we therefore do not interpret
this as a developmental change without depth-matched comparison.

---

## Discussion

These data place all three manganese transporters implicated in human
hypermanganesemia — *SLC39A14* (Tuschl et al., 2016), *SLC30A10* (Tuschl et al.,
2012; Quadri et al., 2012) and *SLC39A8* (Park et al., 2015) — in the zebrafish
RPE, the tissue in which we measure manganese accumulation [CROSS-REF FIG. 2A].
This is consistent with prior work showing that *SLC39A8* and *SLC39A14* are
highly expressed in mammalian RPE and that manganese is required for normal
retinal ultrastructure and function; *slc39a14* null zebrafish show dysregulation
of more than thirty phototransduction genes together with impaired visual
background adaptation and an altered optokinetic response (Tuschl et al., 2022).
The RPE is therefore a plausible site at which transporter loss is converted into
a retinal phenotype, rather than the retina being a passive recipient of systemic
manganese.

[TRANSMISSION ELECTRON MICROSCOPY: STATE WHETHER RPE OR PHOTORECEPTOR
ULTRASTRUCTURAL CHANGES ARE GENOTYPE-, EXPOSURE-, OR INTERACTION-DEPENDENT, WITH
THE EFFECT SIZE AND n PER GROUP. THE TRANSCRIPTOMIC RESULT PREDICTS THAT THE RPE
IS THE AFFECTED COMPARTMENT; SAY EXPLICITLY WHETHER THE EM SUPPORTS, REFINES OR
CONTRADICTS THAT.]

[IMMUNOHISTOCHEMISTRY: STATE WHICH TRANSPORTER PROTEIN WAS LOCALISED, IN WHICH
LAYER, AND WHETHER THE DISTRIBUTION MATCHES THE TRANSCRIPT-LEVEL RESTRICTION
REPORTED HERE. PROTEIN-LEVEL CONFIRMATION IS THE STRONGEST AVAILABLE ANSWER TO
THE DROPOUT LIMITATION BELOW.]

[BULK RNA-SEQ: STATE WHETHER THE TRANSPORTERS AND THE PHOTOTRANSDUCTION MODULE
ARE DIFFERENTIALLY EXPRESSED IN WHOLE EYE, AND RECONCILE DIRECTION AND EFFECT
SIZE WITH THE SINGLE-CELL RESULT. A BULK EFFECT DILUTED ACROSS CELL TYPES IS
EXPECTED IF EXPRESSION IS RPE-RESTRICTED, AND THAT PREDICTION IS TESTABLE HERE.]

Three caveats bound the single-cell result. The RPE population is small (246 and
163 cells), and one 8 dpf library contributes 32% of that stage's RPE cells, so
per-stage statistics should be reported per library rather than pooled. Detection
fraction is bounded below by sequencing depth, so absence of detection is not
evidence of absence — this applies particularly to *slc30a10*, which sits close to
the detection floor throughout. Finally, these data are from wild-type animals and
a public atlas not designed for this question; they establish where the
transporters are expressed, not how that expression responds to manganese
challenge, which [CROSS-REF BULK RNA-SEQ / FIG. X] addresses directly.

---

## Figure legends

**Fig. 2C.** UMAP of the ocular compartment (4,572 cells) with the RPE cluster
coloured and all other ocular cells in grey, shown separately for each stage; RPE
cell number per stage is given in the panel title. Cells from both stages share
one embedding, computed on the combined ocular subset.

**Fig. 2D.** Detection of *slc39a14*, *slc39a8* and *slc30a10* across every
ocular cell type at 5 and 8 dpf, bracketed by controls spanning a range of
expression levels: RPE identity markers (*rpe65a*, *tyrp1b*), abundance-matched
positives (*slc11a2*, *best1*), ambient-RNA sentinels (*rho*, *gnat2*, *elavl3*)
and out-of-tissue negatives (*cryaa*, *hbae1.1*). Dot size, fraction of cells in
which the gene is detected, on an absolute scale to 100%; colour, mean
log-normalised expression among expressing cells only. Detection fraction is
shown rather than mean expression because, for low-abundance transcripts, a mean
across all cells conflates abundance with dropout.

**Fig. S3.** Single-cell reanalysis quality control and controls.
**(A)** Marker dot plot for all ocular sub-clusters, showing the evidence behind
each annotation. **(B)** Cells per cluster at 5 and 8 dpf. **(C)** Detection
fraction for the three transporters against positive and negative control genes,
RPE versus all other ocular clusters. **(D)** Acceptance-criteria table,
auto-evaluated from the analysis outputs. Ambient-RNA correction was applied per
library with CellBender; replication in a second public atlas (Sur et al., 2023)
is not included.
