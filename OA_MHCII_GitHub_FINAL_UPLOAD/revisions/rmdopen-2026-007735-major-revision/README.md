# OA synovium antigen-presentation — RMD Open Major Revision

Manuscript: **Antigen-presentation enrichment in osteoarthritis synovium reflects tissue composition, not uniform cellular activation**

Revision identifier: **RMD Open rmdopen-2026-007735**


## Scope

The package documents Reviewer 1 Major Comments 2, 3, 4 and 6 and Minor Comments 1, 2 and 6. The statistical results were not rerun during package assembly. Existing frozen outputs were copied with SHA-256 tracking. Where a script requires the complete project archive or controlled-access data, this README says so explicitly; the upload copy must not be described as a self-contained end-to-end rerun.

## Directory structure

- `scripts/`: nine archived scripts, including composition correlations, MuSiC fit/LODO diagnostics, attenuation, provenance and read-only audit helpers.
- `config/`: fixed ten-gene index, 11-program definitions and module-freeze metadata.
- `metadata/`: public GEO sample metadata, dataset inventory and the GSE12021 crosswalk.
- `supplementary/`: aggregate reviewer-requested evidence tables.
- `provenance/`: SHA-256 manifest, source mapping and confirmed R environment records.
- `docs/`: evidence report, Supplementary Methods additions, crosswalk, limitations and release checklist.

## Fixed definitions

The fixed antigen-presentation index contains HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, CD74, HLA-DMA, HLA-DMB and CIITA. The separate comparator panel contains 11 prespecified programs; its 60% and ≥3-gene sample-coverage rule does not define a minimum threshold for the primary ten-gene index.

## Reproduction status

- Reviewer 1 Major 2–4 numerical outputs: traceable to frozen result tables and archived scripts; rerunning requires the full project archive and R dependencies.
- Reviewer 1 Major 6 and Minor 1/2/6: audit outputs are available; some platform annotation provenance and patient-level independence remain unresolved.
- MuSiC diagnostics retain native negative `r.squared.full` values. Reference A and B-v2 are OA-derived references.
- Reference A’s 25,021-cell deconvolution branch and the 25,938-cell broad GSE248453 branch are different frozen inputs; the 917-cell cell-level provenance difference is not fully mapped.

## Known limitations

The 158 records are not 158 confirmed independent patients. GSE283079/GSE283080 patient-level independence is unconfirmed. Some original platform annotation provenance is incomplete. The integrated package therefore supports transparent audit and aggregate-result review, not a claim of complete independent end-to-end rerun from this folder alone.

## Software

Confirmed R records: R 4.5.2, metafor 4.8-0 and Matrix 1.7-4. Python package versions were not captured in the archived Major Revision logs and remain UNVERIFIED. See `provenance/environment_versions.md`.

## Manual upload note

The version namespace prevents root-level filename collisions when the complete `revisions/` folder is uploaded. Review `FINAL_CONFLICT_REPORT.md` and `provenance/CONFLICT_CHECK.csv` before uploading.
