# Final Upload Instructions

## Folder to upload

In the GitHub web uploader, drag the complete folder:

`revisions/`

GitHub should then create this path:

`revisions/rmdopen-2026-007735-major-revision/`

Do not drag the contents of `rmdopen-2026-007735-major-revision/` directly into the repository root.

## Expected paths

The uploaded revision should contain `README.md`, `DATA_ACCESS.md`, `scripts/`, `config/`, `metadata/`, `supplementary/`, `provenance/`, and `docs/` beneath the versioned directory. The version directory contains 54 files.

## Verification

After upload, compare the GitHub file tree with the local directory and confirm that the 9 scripts, all Supplementary CSV files, `provenance/CONFLICT_CHECK.csv`, and `provenance/SHA256_UPLOAD_MANIFEST.csv` are visible. The manifest covers 53 files and intentionally excludes itself to avoid circular hashing.

To obtain the complete commit SHA, open the new commit in GitHub and copy the full 40-character SHA from the commit details page; do not infer it from the short web label.

## Error checks

Stop and remove the pending web-upload changes if GitHub shows root-level replacement of `README.md`, `DATA_ACCESS.md`, `scripts/`, or another existing directory, or if an extra `OA_MHCII_GitHub_FINAL_UPLOAD/` directory appears in the repository.

No commit, push, tag, Release, or remote-repository modification was performed by this task.
