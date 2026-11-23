# GitHub release checklist — no push performed

Target repository: https://github.com/Mrxiaohe123/oa-mhcii-apc-synovial-axis

## Staged in this package
- `github_ready_revision/FILE_PROVENANCE_MANIFEST.csv`
- reviewer-analysis scripts under `github_ready_revision/scripts/`
- fixed gene/module definitions and dataset/sample metadata under `github_ready_revision/inputs/`
- selected aggregate reviewer outputs under `github_ready_revision/outputs/`

## Before publication
- Review every staged file and remove any restricted or participant-level data.
- Add all Major Revision scripts and the exact frozen input manifests used by the final manuscript.
- Add R/Python version records and package lock/session files.
- Add README commands distinguishing full reruns from archive-only verification.
- Verify SHA-256 manifest after final staging.
- Commit without rewriting history; obtain the real commit SHA.
- Create a tagged release or DOI only after author approval.
- Update manuscript Data/Code Availability with the real immutable identifier.

No `git push`, release creation, force operation or history rewrite was performed.
