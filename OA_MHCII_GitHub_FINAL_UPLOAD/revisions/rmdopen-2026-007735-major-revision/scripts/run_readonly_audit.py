#!/usr/bin/env python3
"""Read-only reproducibility audit for RMD Open Reviewer 1 minor comments 1, 2 and 6.

This script reads frozen/project-archived files and writes only to the new audit
directory. It performs gene-set presence accounting and file/package inventory;
it does not rerun statistical analyses or modify existing files.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revision_R1_minor1_2_6_audit"
BASE = ROOT / "全部分析内容"
EXPR = BASE / "phase4f_reference_based_deconvolution_sensitivity/03_deconvolution_inputs"
MODULE = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "HLA-DMA", "HLA-DMB", "CIITA"]
DATASETS = ["GSE1919", "GSE12021", "GSE55235", "GSE55457", "GSE82107", "GSE206848", "GSE89408", "GSE143514"]

def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))

def gene_set(path: Path):
    out = set()
    with path.open(encoding="utf-8", errors="replace") as fh:
        header = next(fh, "")
        for line in fh:
            value = line.rstrip("\r\n").split("\t", 1)[0].strip().upper().replace(" ", "")
            if value:
                out.add(value)
    return out

def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def write_csv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main():
    (OUT / "results").mkdir(parents=True, exist_ok=True)
    (OUT / "reports").mkdir(parents=True, exist_ok=True)
    coverage = []
    for ds in DATASETS:
        path = EXPR / f"{ds}_bulk_expression_for_music.tsv"
        present = sorted(gene_set(path) & set(MODULE), key=MODULE.index)
        missing = [g for g in MODULE if g not in present]
        coverage.append({
            "dataset": ds,
            "requested_gene_count": 10,
            "available_gene_count": len(present),
            "available_genes": ";".join(present),
            "missing_genes": ";".join(missing),
            "primary_score_status": "included_in_archived_gate1_input",
            "main_index_minimum_threshold_in_code": "NO_UNIFIED_THRESHOLD_FOUND",
            "rule_source": "全部分析内容/gate1_composition_linkage/scripts/run_gate1_composition_linkage.R:38-41; 76-92",
            "input_source": str(path.relative_to(ROOT)),
            "input_header_line": 1,
            "note": "Coverage is gene-row availability in the frozen expression input; no comparator 60% rule applied."
        })
    write_csv(OUT / "results/primary_index_coverage_audit.csv", coverage, list(coverage[0]))

    inv = {r["accession_id"]: r for r in read_csv(BASE / "manuscript_master/10_bioinformatics_addendum/v2.1_bioinformatics_addendum/01_dataset_flow_inventory/Supplementary_Table_S1_Dataset_Inventory_v2.1.csv")}
    tech_meta = {
        "GSE1919": ("GPL91", "prior_processed_gene_level", "gene-level project matrix; archived input is used as supplied", "No probe mapping in current bridge script; prior processed input", "duplicate symbols averaged after log2(x+1) in the primary bridge script", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:354-365"),
        "GSE12021": ("GPL96 canonical; GPL97 excluded", "series_matrix_plus_platform_annotation", "probe intensity -> mean per mapped gene -> log2(mean+1); GPL96 only in corrected primary", "GPL96.csv/GPL97.csv; first symbol before ///,//,; or comma; uppercase/trim", "mapped probes averaged per gene; corrected script also requires one GPL96 profile per patient", "全部分析内容/gse12021_independent_unit_correction/scripts/run_GSE12021_correction.R:111-145,165-184; :75-109"),
        "GSE55235": ("GPL96", "prior_processed_gene_level", "gene-level project matrix; log2(x+1) in bridge preparation", "UNVERIFIED in current audit beyond prior processed input", "duplicate symbols averaged after log2(x+1)", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:354-365,505-507"),
        "GSE55457": ("GPL96", "prior_processed_gene_level", "gene-level project matrix; log2(x+1) in bridge preparation", "UNVERIFIED in current audit beyond prior processed input", "duplicate symbols averaged after log2(x+1)", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:354-365,505-507"),
        "GSE82107": ("GPL570", "series_matrix_plus_platform_annotation", "series matrix numeric values -> gene mapping -> log2(mean+1)", "GPL570.csv; selected annotation columns are parsed by the shared function", "mapped probes averaged per gene after mapping", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:395-412; :319-340"),
        "GSE206848": ("GPL570", "supplementary_dchip_expression", "dChip expression values -> gene mapping -> log2(mean+1)", "GPL570.csv", "mapped probes averaged per gene after mapping", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:415-455; :319-340"),
        "GSE89408": ("GPL11154", "supplementary_count_matrix", "count matrix -> log2(count+1); no cross-cohort pooled normalization", "UNVERIFIED: gene symbols are read directly from the count matrix; no platform annotation call in this function", "duplicate gene symbols averaged after log2(count+1)", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:458-472"),
        "GSE143514": ("GPL20795", "supplementary_mrna_raw_count", "raw-count matrix -> log2(count+1); GPL20795 records selected", "UNVERIFIED: gene symbols are read directly from count matrix; no probe annotation call in this function", "duplicate gene symbols averaged after log2(count+1)", "全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:475-502"),
    }
    tech = []
    for ds in DATASETS:
        plat, source, scale, annot, dup, code = tech_meta[ds]
        row = inv.get(ds, {})
        # The v2.1 inventory retains the historical 20/13 GSE12021 platform
        # record counts; the canonical corrected analysis is 10/9 GPL96
        # independent profiles, as documented in the correction package.
        oa_n = "10" if ds == "GSE12021" else row.get("OA_n", "UNVERIFIED")
        ctrl_n = "9" if ds == "GSE12021" else row.get("control_n", "UNVERIFIED")
        tech.append({
            "dataset": ds, "platform": plat, "canonical_OA_n": oa_n, "canonical_control_n": ctrl_n,
            "expression_input": f"全部分析内容/phase4f_reference_based_deconvolution_sensitivity/03_deconvolution_inputs/{ds}_bulk_expression_for_music.tsv",
            "source_type": source, "expression_scale_and_transform": scale, "annotation_and_probe_mapping": annot,
            "duplicate_gene_handling": dup, "index_available_gene_count": next(x["available_gene_count"] for x in coverage if x["dataset"] == ds),
            "index_missing_genes": next(x["missing_genes"] for x in coverage if x["dataset"] == ds), "evidence_file_and_code_lines": code,
            "status": "PASS_WITH_UNVERIFIED_FIELDS" if "UNVERIFIED" in " ".join([annot, scale]) else "PASS"
        })
    write_csv(OUT / "results/supplementary_technical_table.csv", tech, list(tech[0]))

    public_root = BASE / "OA_MHCII_APC_GitHub_and_Submission_Package_v1/02_GITHUB_REPOSITORY_READY"
    checks = [
        ("public GitHub URL or fixed remote", "UNVERIFIED", "No URL/remote found in the local archive; workspace is not a Git repository."),
        ("public-ready bulk script", "PRESENT_IN_LOCAL_STAGING_ONLY", "02_GITHUB_REPOSITORY_READY/03_scripts/01_bulk_disease_state_bridge/01_01_bulk_disease_state_bridge.py"),
        ("GSE12021 correction script", "MISSING_FROM_PUBLIC_READY_BUNDLE", "Present in 05_REPRODUCIBILITY/CODE_RELEASE and source archive; not listed in public-ready script inventory."),
        ("frozen primary input manifests/source tables", "PARTIAL", "Public-ready bundle contains selected outputs/tables, but exact corrected GSE12021 crosswalk and all phase4f input manifests are not verified as public."),
        ("Major Revision outputs", "NOT_VERIFIED_PUBLIC", "revision_R1_major2/3/4/6 are outside the staged public-ready tree."),
        ("figure source data", "PARTIAL", "Selected figures and output inventories are staged; complete source-driven Figure 1-6 evidence cannot be confirmed from a public URL."),
        ("environment/session records", "PARTIAL", "Local session logs exist; public-ready package has no confirmed fixed environment lock/commit or public URL."),
        ("release/commit identifier", "UNVERIFIED", "No actual public repository URL, commit, tag, or release was confirmed; placeholders remain in data/code availability files."),
    ]
    write_csv(OUT / "results/public_release_audit.csv", [{"item":a,"status":b,"evidence":c} for a,b,c in checks], ["item","status","evidence"])

    report = f'''# RMD Open Reviewer 1 Minor Comments 1, 2 and 6 — Read-only reproducibility audit\n\n## Scope and status\n\nThis audit reads archived scripts, manifests, frozen input tables and the staged public-ready bundle. No statistical model was rerun and no frozen file was edited. New outputs are confined to `{OUT.relative_to(ROOT)}/`.\n\n## Confirmed facts\n\n- The fixed primary index contains 10 genes: {', '.join(MODULE)}.\n- The Gate 1 scoring function uses the intersection of requested genes and expression rows, then averages gene-wise standardized values; the inspected code does not define a unified minimum coverage threshold for the primary index (`gate1...run_gate1_composition_linkage.R:38-41`).\n- Archived observed coverage is 8/10 for GSE1919 (missing HLA-DRB1 and HLA-DQA1) and 10/10 for the other seven primary datasets.\n- The 60% and >=3 rule belongs to the Gate 4 eleven-program comparator eligibility workflow, not the primary Gate 1 index.\n- GSE12021 correction is explicit: one GPL96 profile per patient, 10 OA/9 control; GPL97 has 1/10 index genes and is technically ineligible.\n\n## Verdicts\n\n| Audit | Verdict | Reason |\n|---|---|---|\n| Primary ten-gene coverage rule | MISMATCH | Manuscript-level wording may imply a prespecified coverage rule, but the inspected Gate 1 scoring code has no unified primary minimum; observed coverage is heterogeneous. |\n| Platform preprocessing | UNVERIFIED | Several dataset branches are traceable, but annotation provenance for prior processed matrices and direct gene-symbol count matrices is incomplete in the inspected archive. GSE12021 is the strongest fully traceable branch. |\n| Code/version/public release | UNVERIFIED | A GitHub-ready staging directory exists, but no public URL, remote, fixed commit, tag or release was confirmed; Major Revision outputs are not verified in that staged bundle. |\n\n## Methods-ready shortest sentences\n\n1. “For the primary ten-gene index, we used the genes present in each frozen dataset-level expression matrix and calculated an equal-weight mean of gene-wise standardized values; no additional universal minimum-coverage threshold was identified in the archived primary scoring code, and dataset-specific coverage was recorded.”\n2. “Platform-specific preprocessing followed the frozen input branches: probe-level matrices were mapped through the corresponding platform annotation and duplicate mapped genes were averaged, whereas gene-level count or processed matrices were transformed as specified in the source workflow; no cross-cohort pooling of raw expression values was performed.”\n3. “GSE12021 was restricted to one GPL96 profile per GEO patient (10 osteoarthritis and 9 control participants); GPL97 records were excluded from the primary index because only 1 of 10 index genes was available.”\n4. “Analysis scripts and derived outputs are retained in the project archive; a public repository URL and immutable release identifier were not available at the time of this audit.”\n\n## Minimal public release checklist\n\nBefore claiming public reproducibility, release or archive: the exact GSE12021 correction script and crosswalk; all eight frozen bulk expression-input manifests and metadata; frozen ten-gene definition and hashes; Gate 1/2/3/4 and Major Revision scripts used in the final manuscript; source-driven Figure 1–6 tables; Supplementary tables; session/package lockfiles; README with relative-path commands; and a DOI or fixed Git tag/commit.\n'''
    (OUT / "reports/minor_comments_1_2_6_audit.md").write_text(report, encoding="utf-8")
    evidence = '''# Evidence map\n\n## Minor Comment 1 — primary ten-gene index coverage\n\n- Gene definition: `全部分析内容/manuscript_master/07_supplement/Supplementary_Table_Module_Freeze_Metadata_v1.5.1.csv`, first data row; `全部分析内容/phase4a_scrna_myeloid_localization/00_governance/myis_oa_tx_frozen_geneset.yaml`.\n- Scoring logic: `全部分析内容/gate1_composition_linkage/scripts/run_gate1_composition_linkage.R:38-41, 76-92`.\n- Dataset-level observed coverage: `results/primary_index_coverage_audit.csv`; archived Gate 1 corroboration: `全部分析内容/gate1_composition_linkage/eligibility_table.csv`.\n- Comparator rule: `全部分析内容/gate4_program_composition_specificity/scripts/validate_gate4.R:23-24, 36-38` and `audit_gate4_inputs.R:7`.\n\n## Minor Comment 2 — platform preprocessing\n\n- Shared preprocessing and platform branches: `全部分析内容/layer25s_disease_state_recovery_and_symptom_boundary/run_layer25s_disease_state_recovery_and_symptom_boundary.py:319-340, 354-365, 368-412, 415-455, 458-502, 505-518`.\n- GSE12021 correction: `全部分析内容/gse12021_independent_unit_correction/scripts/run_GSE12021_correction.R:75-145, 165-184`; coverage tables `GSE12021_GPL96_module_coverage.csv` and `GSE12021_GPL97_coverage_and_sensitivity.csv`; README `README_GSE12021_correction.md`.\n- Supplementary table: `results/supplementary_technical_table.csv`.\n\n## Minor Comment 6 — code/version/public availability\n\n- Local code release manifest: `05_REPRODUCIBILITY/CODE_RELEASE/CODE_RELEASE_MANIFEST.csv`.\n- Local release description: `05_REPRODUCIBILITY/CODE_RELEASE/CODE_RELEASE_README.md`.\n- GitHub-ready staging: `全部分析内容/OA_MHCII_APC_GitHub_and_Submission_Package_v1/02_GITHUB_REPOSITORY_READY/README.md`.\n- Release readiness: `.../09_reproducibility_reports/github_release_readiness_report.md`.\n- Explicit pending URL/DOI placeholders: `.../01_JOURNAL_SUBMISSION_FILES/09_Data_Code_Availability/09_Data_Code_Availability_OA_MHCII_APC.md`.\n- No workspace Git remote: `git remote -v` returned “not a git repository”.\n'''
    (OUT / "reports/minor_comment_evidence_map.md").write_text(evidence, encoding="utf-8")
    (OUT / "reports/methods_ready_sentences.md").write_text(report.split("## Methods-ready shortest sentences",1)[1].split("## Minimal public release checklist",1)[0].strip()+"\n", encoding="utf-8")

if __name__ == "__main__":
    main()
