#!/usr/bin/env python3
"""Assemble an evidence-integration audit from existing frozen outputs.

Read-only with respect to the existing project: this script only copies selected
non-sensitive aggregate files/scripts into revision_evidence_integration/ and
writes audit tables and reports there. It does not execute analysis models.
"""
from __future__ import annotations
import csv, hashlib, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
OUT = ROOT

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def copy(src: str, dest: str, rows, role: str, status: str = 'VERIFIED'):
    p = PROJECT / src
    q = OUT / dest
    q.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        shutil.copy2(p, q)
        rows.append({'source_file': src, 'staged_file': dest, 'sha256': sha(p), 'status': status, 'role': role})
    else:
        rows.append({'source_file': src, 'staged_file': dest, 'sha256': '', 'status': 'UNVERIFIED_MISSING_SOURCE', 'role': role})

def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ['status']
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

def main():
    files = []
    # Existing reviewer analyses and the prior minor-comment audit.
    selected = [
      ('revision_R1_major2/results/cohort_level_correlations.csv','supplementary/Major2_cohort_level_correlations.csv','R1 Major 2'),
      ('revision_R1_major2/results/group_specific_meta_results.csv','supplementary/Major2_group_specific_meta_results.csv','R1 Major 2'),
      ('revision_R1_major2/results/sample_matching_audit.csv','supplementary/Major2_sample_matching_audit.csv','R1 Major 2'),
      ('revision_R1_major3/final_audit/Supplementary_Table_Major3_fit_quality.csv','supplementary/Major3_MuSiC_fit_quality.csv','R1 Major 3'),
      ('revision_R1_major3/final_audit/Supplementary_Table_Major3_LODO_correlations.csv','supplementary/Major3_LODO_correlations.csv','R1 Major 3'),
      ('revision_R1_major3/final_audit/Supplementary_Table_Major3_LODO_figure3.csv','supplementary/Major3_LODO_figure3.csv','R1 Major 3'),
      ('revision_R1_major4/results/mean_attenuation_confidence_intervals.csv','supplementary/Major4_mean_attenuation_CI.csv','R1 Major 4'),
      ('revision_R1_major4/results/cohort_level_attenuation.csv','supplementary/Major4_cohort_attenuation.csv','R1 Major 4'),
      ('revision_R1_major4/results/collinearity_diagnostics.csv','supplementary/Major4_collinearity_diagnostics.csv','R1 Major 4'),
      ('revision_R1_major4/results/LODO_attenuation_stability.csv','supplementary/Major4_LODO_attenuation.csv','R1 Major 4'),
      ('revision_R1_major6/results/dataset_provenance_audit.csv','supplementary/Major6_dataset_provenance.csv','R1 Major 6'),
      ('revision_R1_major6/results/sample_independence_matrix.csv','supplementary/Major6_sample_independence.csv','R1 Major 6'),
      ('revision_R1_major6/results/holdout_reference_independence.csv','supplementary/Major6_holdout_reference_independence.csv','R1 Major 6'),
      ('revision_R1_major6/results/cross_dataset_overlap_candidates.csv','supplementary/Major6_overlap_candidates.csv','R1 Major 6'),
      ('revision_R1_minor1_2_6_audit/results/primary_index_coverage_audit.csv','supplementary/Primary_index_coverage.csv','R1 Minor 1'),
      ('revision_R1_minor1_2_6_audit/results/supplementary_technical_table.csv','supplementary/Platform_technical_table.csv','R1 Minor 2'),
      ('revision_R1_minor1_2_6_audit/results/public_release_audit.csv','supplementary/Public_release_audit.csv','R1 Minor 6'),
    ]
    for a,b,c in selected: copy(a,b,files,c)

    # Read-only staging of real scripts/configs required for reviewer analyses.
    scripts = [
      ('全部分析内容/gse12021_independent_unit_correction/scripts/run_GSE12021_correction.R','scripts/run_GSE12021_correction.R','canonical GSE12021 correction'),
      ('全部分析内容/gate1_composition_linkage/scripts/run_gate1_composition_linkage.R','scripts/run_gate1_composition_linkage.R','Gate 1 primary index/composition'),
      ('revision_R1_major2/scripts/run_major_comment2.R','scripts/run_major_comment2.R','R1 Major 2'),
      ('revision_R1_major3/scripts/run_major_comment3_LODO.R','scripts/run_major_comment3_LODO.R','R1 Major 3'),
      ('revision_R1_major4/scripts/run_major_comment4_audit.R','scripts/run_major_comment4_audit.R','R1 Major 4'),
      ('revision_R1_major6/scripts/run_major_comment6_provenance_audit.py','scripts/run_major_comment6_provenance_audit.py','R1 Major 6'),
      ('revision_R1_minor1_2_6_audit/scripts/run_readonly_audit.py','scripts/run_readonly_audit.py','R1 Minor 1/2/6'),
      ('全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py','scripts/run_layer25s_bridge.py','primary preprocessing branches'),
    ]
    for a,b,c in scripts: copy(a,b,files,c)
    configs = [
      ('全部分析内容/phase4a_scrna_myeloid_localization/00_governance/myis_oa_tx_frozen_geneset.yaml','inputs/frozen_primary_gene_set.yaml','fixed module definition'),
      ('全部分析内容/manuscript_master/07_supplement/Supplementary_Table_Module_Freeze_Metadata_v1.5.1.csv','inputs/module_freeze_metadata.csv','module provenance'),
      ('全部分析内容/gate4_program_composition_specificity/program_definitions.csv','inputs/program_definitions.csv','11-program definitions'),
      ('全部分析内容/manuscript_master/10_bioinformatics_addendum/v2.1_bioinformatics_addendum/01_dataset_flow_inventory/Supplementary_Table_S1_Dataset_Inventory_v2.1.csv','inputs/dataset_inventory.csv','dataset inventory'),
      ('全部分析内容/gse12021_independent_unit_correction/source_data/GSE12021_source_record_crosswalk.csv','inputs/GSE12021_source_record_crosswalk.csv','GSE12021 sample crosswalk'),
    ]
    for a,b,c in configs: copy(a,b,files,c)
    for p in sorted((PROJECT/'全部分析内容/phase4f_reference_based_deconvolution_sensitivity/03_deconvolution_inputs').glob('*_metadata_for_deconvolution.csv')):
        copy(str(p.relative_to(PROJECT)), 'inputs/'+p.name, files, 'canonical sample metadata')

    write_csv(OUT/'FILE_PROVENANCE_MANIFEST.csv', files)
    claim_rows = [
      {'Reviewer_comment':'R1 Minor 1','Manuscript_section':'Methods / fixed index','Exact_claim_or_result':'Primary index uses available genes; GSE1919 8/10 and other seven datasets 10/10; no universal primary threshold in code.','Analysis_script':'全部分析内容/gate1_composition_linkage/scripts/run_gate1_composition_linkage.R:38-41,76-92','Frozen_input':'phase4f.../03_deconvolution_inputs/*_bulk_expression_for_music.tsv','Derived_output':'revision_R1_minor1_2_6_audit/results/primary_index_coverage_audit.csv','Supplementary_evidence':'supplementary/Primary_index_coverage.csv','Reproduction_status':'VERIFIED','Evidence_limitations':'The manuscript phrase “prespecified coverage rule” is broader than the primary scoring code.','Required_action':'Clarify primary-index rule separately from comparator rule.'},
      {'Reviewer_comment':'R1 Minor 2','Manuscript_section':'Methods / preprocessing','Exact_claim_or_result':'Platform-specific input handling and GSE12021 GPL96 correction.','Analysis_script':'layer25s.../run_layer25s_disease_state_recovery_and_symptom_boundary.py:319-340,354-518; gse12021.../run_GSE12021_correction.R:75-145','Frozen_input':'platform annotation files; GSE12021 crosswalk and coverage tables','Derived_output':'revision_R1_minor1_2_6_audit/results/supplementary_technical_table.csv','Supplementary_evidence':'supplementary/Platform_technical_table.csv','Reproduction_status':'PARTIALLY_VERIFIED','Evidence_limitations':'Prior processed matrices and direct gene-symbol count-matrix provenance remain incomplete for some datasets.','Required_action':'Mark untraceable annotation fields UNVERIFIED or provide source manifests.'},
      {'Reviewer_comment':'R1 Minor 6','Manuscript_section':'Data/code availability','Exact_claim_or_result':'Major Revision scripts and evidence outputs are publicly released and versioned.','Analysis_script':'05_REPRODUCIBILITY/CODE_RELEASE/CODE_RELEASE_MANIFEST.csv; local GitHub-ready staging','Frozen_input':'public-ready package and local archive','Derived_output':'FILE_PROVENANCE_MANIFEST.csv; supplementary/Public_release_audit.csv','Supplementary_evidence':'GitHub-ready staging directory; target repository README','Reproduction_status':'PARTIALLY_VERIFIED','Evidence_limitations':'Public repository is confirmed, but its visible one-commit contents do not yet demonstrate synchronized Major Revision scripts/results; immutable commit was not recovered by local remote query.','Required_action':'Author must review staging, commit, and publish a fixed release before claiming complete public reproducibility.'},
      {'Reviewer_comment':'R1 Major 2','Manuscript_section':'Results / composition correlations','Exact_claim_or_result':'OA-only, control-only and OA-adjusted correlations with REML/HK on Fisher-z scale.','Analysis_script':'revision_R1_major2/scripts/run_major_comment2.R','Frozen_input':'revision_R1_major2/results/canonical_matched_input.csv','Derived_output':'revision_R1_major2/results/group_specific_meta_results.csv','Supplementary_evidence':'supplementary/Major2_group_specific_meta_results.csv','Reproduction_status':'VERIFIED','Evidence_limitations':'Group-specific meta-analysis uses small k and prediction intervals may be unavailable/degenerate for Reference B OA.','Required_action':'Retain exploratory wording and full cohort-level table.'},
      {'Reviewer_comment':'R1 Major 3','Manuscript_section':'Results / deconvolution robustness','Exact_claim_or_result':'Negative MuSiC R2 values retained; 8 Reference A and 4 Reference B-v2 donor LODO runs.','Analysis_script':'revision_R1_major3/final_audit/scripts/run_final_audit.R','Frozen_input':'revision_R1_major3/results/reference_audit.csv and LODO outputs','Derived_output':'revision_R1_major3/final_audit/*.csv','Supplementary_evidence':'supplementary/Major3_MuSiC_fit_quality.csv; Major3_LODO_correlations.csv','Reproduction_status':'VERIFIED','Evidence_limitations':'R2 is a model-fit diagnostic, not direct cell-fraction accuracy.','Required_action':'Do not truncate negative R2 or overstate donor robustness.'},
      {'Reviewer_comment':'R1 Major 4','Manuscript_section':'Results / composition sensitivity','Exact_claim_or_result':'Mean paired attenuation and t/Bootstrap CIs; no causal mediation interpretation.','Analysis_script':'revision_R1_major4/scripts/run_major_comment4_audit.R','Frozen_input':'revision_R1_major4/results/cohort_level_attenuation.csv','Derived_output':'revision_R1_major4/results/mean_attenuation_confidence_intervals.csv','Supplementary_evidence':'supplementary/Major4_mean_attenuation_CI.csv; Major4_collinearity_diagnostics.csv','Reproduction_status':'VERIFIED','Evidence_limitations':'Bootstrap is cohort-resampling only; it does not propagate MuSiC uncertainty.','Required_action':'Report t-based CIs including those crossing zero.'},
      {'Reviewer_comment':'R1 Major 5/7/8 + R2','Manuscript_section':'Discussion / interpretation','Exact_claim_or_result':'Small donor numbers, descriptive 11-program ranking, exploratory clinical context, and coexisting composition/intrinsic activity.','Analysis_script':'Existing frozen manuscript/reviewer audit files; no new analysis requested','Frozen_input':'FINAL_MANUSCRIPT_MASTER.md and frozen Gate3c/Gate4 outputs','Derived_output':'Existing revision reports','Supplementary_evidence':'Crosswalk and methods additions in this integration package','Reproduction_status':'PARTIALLY_VERIFIED','Evidence_limitations':'This integration audit does not rewrite the manuscript or perform new experiments.','Required_action':'Author review of minimal patches and discussion wording.'},
    ]
    write_csv(OUT/'CLAIM_EVIDENCE_MATRIX.csv', claim_rows)

    (OUT/'SUPPLEMENTARY_METHODS_ADDITIONS.md').write_text('''# Supplementary Methods additions — author-review draft\n\n## Primary bulk index and coverage\nThe fixed antigen-presentation index comprised HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, CD74, HLA-DMA, HLA-DMB and CIITA. Within each frozen dataset-level matrix, available genes were standardized across samples and averaged with equal gene weights; missing genes were not imputed. The archived primary scoring code did not define a universal minimum-coverage threshold. Observed coverage was 8/10 in GSE1919 and 10/10 in the other seven primary datasets. The separate 60% and ≥3-gene rule applied only to the prespecified 11-program comparator panel.\n\n## GSE12021 correction\nGSE12021 was restricted to one GPL96 profile per GEO patient identifier (10 OA and 9 controls). GPL97 records were excluded from the canonical primary analysis because only 1/10 fixed-index genes was available. Probe annotations used the locked first-symbol rule; probes mapping to the same gene were averaged on the expression scale and transformed as log2(mean + 1).\n\n## Composition correlations and adjustment\nFor Reviewer 1 Major Comment 2, correlations were calculated at the dataset level within OA, within controls, and after adjustment for OA status using rank transformation and residualization. Fisher-z REML/Hartung–Knapp pooling was applied where prespecified eligibility permitted. These analyses are exploratory composition-sensitivity analyses and do not establish causal mediation or independent cell-intrinsic activation.\n\n## Deconvolution diagnostics\nMuSiC `r.squared.full` values were retained on their native scale, including negative values. Reference A used eight OA donors and 25,021 cells in the deconvolution branch; the broader GSE248453 input contained 25,938 cells. These are different frozen branches, and a complete cell-level provenance map explaining the 917-cell difference was not available.\n\n## Provenance and public release\nThe eight primary cohorts comprise 158 canonical expression-analysis records, not 158 verified independent donors. No exact cross-accession GSM/BioSample duplicates were identified in the audited metadata, but patient-level overlap could not be excluded for most GEO series. GSE283079 and the GSE283080-derived Reference B-v2 belong to the same study family; patient-level independence was not confirmed.\n''', encoding='utf-8')

    (OUT/'SUPPLEMENTARY_CROSSWALK.md').write_text('''# Supplementary crosswalk and numbering proposal\n\nExisting Major 2–6 audit tables should be reused where possible. Proposed labels are provisional and must be reconciled with the current RMD Open supplement numbering.\n\n| Proposed item | Existing source | Status |\n|---|---|---|\n| Table R1-2A | Major 2 cohort-level correlations and group-specific meta-analysis | ready for author review |\n| Table R1-3A | Major 3 MuSiC fit quality including negative R² | ready for author review |\n| Table R1-3B | Major 3 donor LODO correlations and Figure 3 robustness | ready for author review |\n| Table R1-4A | Major 4 paired attenuation, t-based and Bootstrap CIs | ready for author review |\n| Table R1-4B | Major 4 collinearity and leave-one-cohort-out diagnostics | ready for author review |\n| Table R1-6A | Provenance, overlap and holdout/reference independence | ready for author review |\n| Table R1-M1/M2 | Primary index coverage and platform preprocessing | ready for author review |\n\nNo existing Supplementary number was overwritten or assigned as final in this staging package.\n''', encoding='utf-8')

    (OUT/'MANUSCRIPT_MINIMAL_PATCHES.md').write_text('''# Manuscript minimal patches — author review only\n\nNo manuscript file was modified.\n\n1. Replace any statement implying that 158 records are 158 independent donors with “158 canonical expression-analysis records; donor-level independence was directly resolved only for GSE12021.”\n2. In the fixed-index Methods paragraph, separate the available-gene scoring rule from the comparator-panel 60%/≥3-gene eligibility rule.\n3. Describe Reference A and Reference B-v2 associations as composition-sensitivity analyses; do not describe attenuation as mediation or as a quantified cellular contribution.\n4. State that negative MuSiC R² values were retained and that R² is a model-fit diagnostic, not direct validation of inferred fractions.\n5. State that GSE283079 was an external holdout for the bulk analysis, while patient-level independence from GSE283080-derived Reference B-v2 was not confirmed.\n6. Preserve cautious single-cell language: “no conclusive evidence of resident-family-specific enrichment” and “insufficient donor-level power to exclude state-specific differences.”\n\nEach patch requires author confirmation against the current RMD Open manuscript copy before insertion.\n''', encoding='utf-8')

    (OUT/'UNVERIFIED_ITEMS.md').write_text('''# Unverified or partially verified items\n\n- No complete cell-level provenance map explains the 917-cell difference between 25,938 broad GSE248453 cells and 25,021 Reference A deconvolution cells.\n- For several prior-processed or direct gene-symbol count matrices, the original annotation/probe-level provenance is not fully recoverable from the current archive.\n- Patient-level overlap across most GEO accessions remains unresolved despite no exact GSM/BioSample duplicate in available metadata.\n- GSE283079 and GSE283080 are in the same study family; donor-level independence is not confirmed.\n- The public GitHub repository exists and is public, but the current visible repository state is not demonstrated to contain the complete Major Revision evidence package.\n- No immutable commit/tag/release for the integrated revision package was confirmed in this audit.\n- Supplementary numbering remains provisional.\n''', encoding='utf-8')

    (OUT/'GITHUB_RELEASE_CHECKLIST.md').write_text('''# GitHub release checklist — no push performed\n\nTarget repository: https://github.com/Mrxiaohe123/oa-mhcii-apc-synovial-axis\n\n## Staged in this package\n- `github_ready_revision/FILE_PROVENANCE_MANIFEST.csv`\n- reviewer-analysis scripts under `github_ready_revision/scripts/`\n- fixed gene/module definitions and dataset/sample metadata under `github_ready_revision/inputs/`\n- selected aggregate reviewer outputs under `github_ready_revision/outputs/`\n\n## Before publication\n- Review every staged file and remove any restricted or participant-level data.\n- Add all Major Revision scripts and the exact frozen input manifests used by the final manuscript.\n- Add R/Python version records and package lock/session files.\n- Add README commands distinguishing full reruns from archive-only verification.\n- Verify SHA-256 manifest after final staging.\n- Commit without rewriting history; obtain the real commit SHA.\n- Create a tagged release or DOI only after author approval.\n- Update manuscript Data/Code Availability with the real immutable identifier.\n\nNo `git push`, release creation, force operation or history rewrite was performed.\n''', encoding='utf-8')

    readme = '''# GitHub-ready revision staging (author review only)\n\nThis directory is a non-pushed staging bundle assembled from existing project files. It contains selected scripts, non-sensitive definitions/metadata, aggregate reviewer outputs and SHA-256 provenance. It does not contain OAI/NDA participant-level data or credentials.\n\nThe bundle is not a public release and has no commit/tag/DOI assigned. Review `FILE_PROVENANCE_MANIFEST.csv` and `../GITHUB_RELEASE_CHECKLIST.md` before publication.\n'''
    (OUT/'github_ready_revision').mkdir(parents=True, exist_ok=True)
    (OUT/'github_ready_revision/README.md').write_text(readme, encoding='utf-8')
    for group in ('scripts', 'inputs', 'supplementary'):
        src_dir = OUT / group
        dst_dir = OUT / 'github_ready_revision' / group
        dst_dir.mkdir(parents=True, exist_ok=True)
        for src in src_dir.iterdir():
            if src.is_file():
                shutil.copy2(src, dst_dir / src.name)

if __name__ == '__main__': main()
