# RMD Open Major Revision — Evidence Integration Report

## Scope

This is an author-review evidence audit and staging package. It reused existing Reviewer 1 Major 2–6 outputs and the prior Minor 1/2/6 audit. No statistical model was rerun, no analysis method was changed, and no frozen manuscript, input, figure or result was overwritten.

## Confirmed evidence

1. The eight primary bulk cohorts contain 158 canonical expression-analysis records, not 158 verified independent donors. GSE12021 is the only cohort with an archived patient-level correction: 19 GPL96 profiles, 10 OA and 9 controls; 14 GPL97 repeat-platform records were excluded.
2. The fixed ten-gene index is unchanged. The archived Gate 1 score function uses available genes and equal-weight gene-wise standardized values. No universal primary-index minimum coverage threshold is implemented in the inspected code. GSE1919 has 8/10 genes, missing HLA-DRB1 and HLA-DQA1; the other seven primary datasets have 10/10. The 60% plus ≥3-gene rule is comparator-only.
3. Reviewer 1 Major 2 outputs support positive pooled OA-only and OA-adjusted correlations, while control-only estimates are less precise and prediction intervals are wide. Full cohort-level results are staged in `supplementary/`.
4. Reviewer 1 Major 3 outputs retain native negative MuSiC `r.squared.full` values and include all 8 Reference A and 4 Reference B-v2 donor-omission runs. These are model-fit diagnostics, not direct validation of cell-fraction accuracy.
5. Reviewer 1 Major 4 outputs use paired cohort-level `abs(unadjusted beta)-abs(adjusted beta)` values. The t-based CIs are preserved, including intervals crossing zero. Bootstrap was cohort-resampling only and does not propagate deconvolution measurement uncertainty.
6. Reviewer 1 Major 6 found no exact cross-accession GSM/BioSample duplicates in the audited metadata, but patient-level overlap remains unresolved for most GEO series. GSE283079 and the GSE283080-derived Reference B-v2 are in the same study family and patient-level independence was not confirmed.
7. Reference A’s 25,021-cell deconvolution input and the 25,938-cell broad GSE248453 input are different frozen branches. The complete 917-cell cell-level provenance map is not available; no specific exclusion explanation was invented.

## Evidence status

The detailed claim-to-source mapping is in `CLAIM_EVIDENCE_MATRIX.csv`. The principal statuses are:

- Major 2, Major 3 and Major 4 numerical outputs: `VERIFIED` against existing frozen result tables and scripts.
- Minor 1: `VERIFIED` for observed coverage and code behavior; documentation needs clarification because the manuscript’s phrase “prespecified coverage rule” is broader than the primary scoring code.
- Minor 2: `PARTIALLY_VERIFIED`; preprocessing branches are traceable, but some prior-processed and direct count-matrix annotation provenance remains incomplete.
- Minor 6: `PARTIALLY_VERIFIED`; the target public repository is real, but the complete Major Revision package and immutable release identifier were not confirmed.

## Public release status

The target repository is public at <https://github.com/Mrxiaohe123/oa-mhcii-apc-synovial-axis>. The repository page identifies it as public and currently displays one commit, but the visible contents do not demonstrate synchronization of the full Major Revision scripts, manifests and reviewer evidence tables. The local integration bundle is therefore a staging package, not a claimed public release. No push, force operation, history rewrite, tag or release was performed.

## Author decisions required before submission

1. Approve the primary-index wording that explicitly separates available-gene scoring from comparator eligibility.
2. Decide whether the remaining UNVERIFIED platform fields can be supported by additional source manifests.
3. Review the provisional Supplementary numbering and cross-references.
4. Review and commit the staged public files, obtain the real immutable commit/tag/DOI, and update Data/Code Availability only afterward.
5. Confirm that all Major Revision response letters and manuscript claims use composition-sensitivity language rather than mediation or independent-donor language.
