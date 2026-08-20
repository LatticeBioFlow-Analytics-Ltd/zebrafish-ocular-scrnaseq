# Results and Discussion — integrated draft (DMM)

Covers Fig. 2C–D, Fig. S3, and cross-references to the bulk RNA-seq (Fig. 2A–B,
3, 4), EM (Fig. 6) and IHC datasets. Every number traces to a file; see
`PROVENANCE.md`. The single unquantified dataset is the IHC, and it is marked as
such wherever it appears.

---

## Results

### Manganese transporters are restricted to the retinal pigment epithelium

To ask whether the eye expresses the machinery to handle manganese locally, we
reanalysed a public single-cell atlas of the zebrafish head at 5 and 8 dpf
(Raj et al., 2020), recovering 107,444 cells in 34 clusters and sub-clustering
the ocular compartment (4,572 cells; Fig. S3A). Ambient RNA was removed per
library before analysis (Fleming et al., 2023), which is necessary here because
the fraction of reads in cells differs systematically between stages
(0.49–0.56 at 5 dpf; 0.33–0.47 at 8 dpf).

A single population expressing *rpe65a*, *tyrp1b* and *pmela* was identified as
retinal pigment epithelium (RPE) at both stages (246 and 163 cells; Fig. 2C).
All three transporters implicated in human hypermanganesemia were enriched in
this population relative to every other ocular cell type: *slc39a14* in 12.2% of
RPE cells against 3.7% elsewhere, *slc39a8* in 5.4% against 0.07%, and *slc30a10*
in 2.2% against 0.00% (Fig. 2D). *slc39a8* and *slc30a10* were essentially
confined to the RPE. Detection of *slc39a14* matched the abundance-matched
positive control *slc11a2* (12.5%), while out-of-tissue controls remained near
zero in every cluster (*cryaa* ≤3.3%; *hbae1.1* ≤0.6%), indicating the signal is
not residual ambient RNA (Fig. 2D, Fig. S3C). This localisation agrees with the
reported enrichment of SLC39A8 and SLC39A14 in mammalian RPE (Tuschl et al.,
2022).

### Manganese challenge produces a genotype-dependent transcriptional response

Whole-eye RNA-sequencing across a fully crossed genotype × exposure design
(*n* = 3 per group) shows that the transcriptional response to manganese is
almost entirely confined to the mutant. Wild-type animals responded to MnCl₂
with 10 differentially expressed genes in total (4 up, 6 down), and the two
genotypes were near-indistinguishable at baseline (23 genes). Exposure of
*slc39a14*⁻ᐟ⁻ animals changed 1,154 genes (752 up, 402 down), and the genotype ×
exposure interaction term recovered 475 (276 up, 199 down; Fig. 2A–B).

The interaction signature is directional and coherent. Endoplasmic-reticulum
stress and protein folding were increased (GO "response to endoplasmic reticulum
stress" NES = +2.22, *P*<sub>adj</sub> = 6.1 × 10⁻⁷; KEGG "protein processing in
endoplasmic reticulum" NES = +2.33, *P*<sub>adj</sub> = 1.9 × 10⁻¹⁰), led by
*atf3*, *atf5a*, *jun* and *creb5a* (Fig. 3). Phototransduction and photoreceptor
structure were decreased (GO "sensory perception of light stimulus" NES = −2.63,
*P*<sub>adj</sub> = 5.4 × 10⁻¹⁵; "photoreceptor outer segment" NES = −2.30; KEGG
"phototransduction" NES = −2.40, *P*<sub>adj</sub> = 5.0 × 10⁻⁹), led by *pde6h*,
*grk1b* and *lrit1b* (Fig. 4). KEGG "retinol metabolism", the visual cycle
carried out by the RPE, was also decreased (NES = −2.08,
*P*<sub>adj</sub> = 3.0 × 10⁻⁴).

Deconvolution of the same libraries against retinal cell-type signatures
reproduced this at the level of composition: cone and rod photoreceptor
signatures fell in the interaction contrast (log FC −1.47 and −1.35, both
*P*<sub>adj</sub> = 0.009) alongside ciliary marginal zone progenitors (−1.38),
while the microglia/macrophage signature rose (+0.87, *P*<sub>adj</sub> = 0.092).
No cell-type signature changed in exposed wild-type animals.

Notably, none of the three transporters was itself differentially expressed in
any contrast (all *P*<sub>adj</sub> > 0.4; base mean 14.8–55). This is the
expected consequence of RPE-restricted expression: the RPE is 4.3% of the ocular
compartment by cell number, so a transcript confined to it is heavily diluted in
whole-eye RNA. The single-cell and bulk results are therefore consistent rather
than discrepant, and the dilution is itself evidence of restriction.

### Mitochondrial ultrastructure diverges between genotypes under manganese

Blinded hand annotation of photoreceptor inner-segment mitochondria at 25,000×
(the only magnification acquired in all four groups; Fig. 6) shows a baseline
genotype effect and a divergent response to exposure. Cristae density was lower
in *slc39a14*⁻ᐟ⁻ larvae irrespective of exposure (Hedges' *g* = −2.11, 95% CI
−5.02 to −1.07 unexposed; *g* = −1.05, −3.96 to −0.07 exposed). Under manganese,
wild-type mitochondrial area fraction fell (*g* = −1.33, −2.84 to −0.50; 0.195 to
0.061) whereas the mutant's did not (0.116 to 0.200), leaving mutant mitochondria
both larger (mean area *g* = +1.03, 0.25 to 2.27) and occupying more of the field
(*g* = +0.91, 0.22 to 1.84) than wild type under the same challenge.

---

## Discussion

Three datasets converge on the same interpretation: the RPE holds the eye's
manganese-handling capacity, and losing *slc39a14* converts a manganese challenge
that wild-type animals absorb without transcriptional consequence into
photoreceptor stress and dysfunction.

The single-cell data place all three transporters implicated in human
hypermanganesemia — *SLC39A14* (Tuschl et al., 2016), *SLC30A10* (Tuschl et al.,
2012; Quadri et al., 2012) and *SLC39A8* (Park et al., 2015) — in the RPE, and
the bulk data show that eyes lacking *slc39a14* mount an unfolded-protein and
integrated-stress response to manganese while wild-type eyes do essentially
nothing (10 genes). Manganese is a cofactor for ER- and Golgi-resident enzymes
and its dyshomeostasis is a recognised cause of glycosylation and secretory
failure (Park et al., 2015), so an ER-stress signature is the expected proximate
readout of a cell that can no longer buffer the metal.

The damage, however, is not read out in the RPE. The RPE signature is unchanged
in the interaction contrast, whereas cone and rod signatures fall and the
microglial signature rises, and the decreased pathways are dominated by
phototransduction and photoreceptor outer-segment terms. The RPE appears to fail
functionally rather than numerically: "retinol metabolism", the visual cycle it
performs on behalf of photoreceptors, is among the decreased KEGG pathways. This
is consistent with the phenotype previously reported for this mutant, in which
more than thirty phototransduction genes are dysregulated alongside impaired
visual background adaptation and an altered optokinetic response (Tuschl et al.,
2022), and it identifies the RPE as the plausible origin of that phenotype rather
than the retina being a passive recipient of systemic manganese.

The ultrastructural data add a mitochondrial dimension and one caution. The
cristae deficit is present in unexposed mutants, so it is a consequence of
transporter loss rather than of manganese loading — consistent with the
manganese *deficiency* that coexists with hypersensitivity in this model (Tuschl
et al., 2022) and with the manganese dependence of mitochondrial superoxide
dismutase. The exposure effect is better described as a failure to respond than
as accumulation: wild-type animals reduce inner-segment mitochondrial area
fraction threefold under manganese while mutants do not. We emphasise this
because an automated segmentation pilot on the same frames produced the opposite
and much more striking result (area fraction 20.4% in mutant + MnCl₂ against
2.4–8.4% elsewhere, *g* = +2.9), which overlay inspection showed to be an
artefact of detection tracking section staining rather than mitochondrial
content; the hand-annotated values reported here supersede it.

Two limitations bound these conclusions. The EM frames are indexed by group but
not by animal, so *n* refers to fields rather than fish and within-animal
correlation cannot be excluded; the effect sizes should be read as descriptive.
And the single-cell atlas is from unexposed wild-type animals, so it establishes
where the transporters are expressed, not how that expression responds to
challenge. Protein-level localisation of RPE65 and the transporters, and
quantification of cone (GNAT2), microglial (4C4), apoptotic (caspase-3) and
synaptic (Ribeye) markers across the same four groups, are in progress on a
blinded 239-fish immunohistochemistry series (Methods); these will test directly
whether the microglial increase and photoreceptor loss inferred from
deconvolution are visible in tissue.
