# Supplementary Methods additions — author-review draft

## Primary bulk index and coverage
The fixed antigen-presentation index comprised HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, CD74, HLA-DMA, HLA-DMB and CIITA. Within each frozen dataset-level matrix, available genes were standardized across samples and averaged with equal gene weights; missing genes were not imputed. The archived primary scoring code did not define a universal minimum-coverage threshold. Observed coverage was 8/10 in GSE1919 and 10/10 in the other seven primary datasets. The separate 60% and ≥3-gene rule applied only to the prespecified 11-program comparator panel.

## GSE12021 correction
GSE12021 was restricted to one GPL96 profile per GEO patient identifier (10 OA and 9 controls). GPL97 records were excluded from the canonical primary analysis because only 1/10 fixed-index genes was available. Probe annotations used the locked first-symbol rule; probes mapping to the same gene were averaged on the expression scale and transformed as log2(mean + 1).

## Composition correlations and adjustment
For Reviewer 1 Major Comment 2, correlations were calculated at the dataset level within OA, within controls, and after adjustment for OA status using rank transformation and residualization. Fisher-z REML/Hartung–Knapp pooling was applied where prespecified eligibility permitted. These analyses are exploratory composition-sensitivity analyses and do not establish causal mediation or independent cell-intrinsic activation.

## Deconvolution diagnostics
MuSiC `r.squared.full` values were retained on their native scale, including negative values. Reference A used eight OA donors and 25,021 cells in the deconvolution branch; the broader GSE248453 input contained 25,938 cells. These are different frozen branches, and a complete cell-level provenance map explaining the 917-cell difference was not available.

## Provenance and public release
The eight primary cohorts comprise 158 canonical expression-analysis records, not 158 verified independent donors. No exact cross-accession GSM/BioSample duplicates were identified in the audited metadata, but patient-level overlap could not be excluded for most GEO series. GSE283079 and the GSE283080-derived Reference B-v2 belong to the same study family; patient-level independence was not confirmed.
