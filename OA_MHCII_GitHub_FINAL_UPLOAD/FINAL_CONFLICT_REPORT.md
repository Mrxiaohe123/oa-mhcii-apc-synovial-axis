# Final Conflict Report

## Scope

This report audits the versioned upload namespace `revisions/rmdopen-2026-007735-major-revision/` against the visible root-level names of the target repository. No remote repository was modified.

## Resolution

The prior eight root-level conflicts are resolved for this upload method by requiring the author to upload the complete `revisions/` directory as a new top-level directory. The revision files therefore occupy a new path namespace and do not replace `README.md`, `DATA_ACCESS.md`, `scripts/`, or any other existing root-level path. This is path isolation, not evidence that the existing files are identical.

| Prior item | Final status | Verification |
|---|---|---|
| `README.md` | RESOLVED_BY_NAMESPACE | Final path is `revisions/rmdopen-2026-007735-major-revision/README.md`; root `README.md` is not included as a replacement. |
| `DATA_ACCESS.md` | RESOLVED_BY_NAMESPACE | Final path is versioned; the previously confirmed root conflict is not overwritten. |
| `scripts/` | RESOLVED_BY_NAMESPACE | Final path is `revisions/rmdopen-2026-007735-major-revision/scripts/`. |
| `config/` | RESOLVED_BY_NAMESPACE | Final path is versioned. |
| `metadata/` | RESOLVED_BY_NAMESPACE | Final path is versioned. |
| `supplementary/` | RESOLVED_BY_NAMESPACE | Final path is versioned. |
| `provenance/` | RESOLVED_BY_NAMESPACE | Final path is versioned. |
| `docs/` | RESOLVED_BY_NAMESPACE | Final path is versioned. |

## Upload condition

Do not drag the contents of the version directory into the repository root. Drag the `revisions/` folder so GitHub creates the complete versioned path.

## Integrity

The version directory contains 54 files. `provenance/SHA256_UPLOAD_MANIFEST.csv` contains hashes for the other 53 files and is excluded from its own hash calculation.
