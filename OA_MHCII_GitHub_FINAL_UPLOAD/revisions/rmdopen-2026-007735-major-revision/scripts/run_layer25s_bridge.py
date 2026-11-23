#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import itertools
import json
import math
import os
import re
import sys
import tarfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "layer25s_disease_state_recovery_and_symptom_boundary"
GOVDIR = OUTDIR / "00_governance"
META1 = OUTDIR / "01_legacy_meta"
HOLD2 = OUTDIR / "02_gse283079_holdout"
LOC3 = OUTDIR / "03_gse216651_localization"
EXP4 = OUTDIR / "04_gse283080_exploratory"
SYN5 = OUTDIR / "05_synthesis"
DATADIR = OUTDIR / "data"
RAWGEO = DATADIR / "raw_geo"
PLATDIR = DATADIR / "platform_annotation"
PROCDIR = DATADIR / "processed"
for path in [OUTDIR, GOVDIR, META1, HOLD2, LOC3, EXP4, SYN5, PROCDIR]:
    path.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
DEPDIR = ROOT / "_pydeps"
if str(DEPDIR) not in sys.path:
    sys.path.insert(0, str(DEPDIR))

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import sparse, stats
from scipy.io import mmread


FROZEN_YAML = ROOT / "layer25r_adult_oa_direct_transcriptomic_bridge" / "01_signature_freeze" / "myis_oa_tx_frozen_geneset.yaml"
PRIOR_BULK = ROOT / "layer25_oa_transcriptomic_bridge_real" / "data" / "bulk"
PRIOR_SYMPTOM = ROOT / "layer25r_v11_audit_and_replication"
PLATFORM_MAP = {
    "GPL96": PLATDIR / "GPL96.csv",
    "GPL97": PLATDIR / "GPL97.csv",
    "GPL570": PLATDIR / "GPL570.csv",
}
LOG_PATH = OUTDIR / "layer25s_run_log.txt"

PRIMARY_MODULE = "MHCII_antigen_presentation_module"
SECONDARY_MODULES = [
    "MyIS_OA_Tx_composite",
    "SPP1_macrophage_module",
    "S100A8_A9_myeloid_module",
    "IL1_IL6_TNF_response_module",
    "complement_coagulation_module",
]
SENSITIVITY_MODULES = ["CNBP_related_module", "IFN_response_module"]
GENERIC_MHCII_REFERENCE = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB",
    "CD74", "CIITA", "B2M", "TAP1", "TAP2", "PSMB8",
]
HOUSEKEEPING_CONTROL = [
    "ACTB", "GAPDH", "HPRT1", "EEF1A1", "PPIA", "RPLP0", "RPL13A", "TBP", "UBC", "SF3B1",
]
CELLTYPE_MARKERS = {
    "myeloid": ["LST1", "FCER1G", "TYROBP", "CTSS", "AIF1", "SAT1"],
    "dendritic": ["FCER1A", "CLEC10A", "CD1C", "HLA-DRA", "HLA-DPA1", "CST3"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "THY1"],
    "lymphocyte": ["CD3D", "TRAC", "NKG7", "MS4A1", "CD79A", "LTB"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "PLVAP"],
    "mast": ["TPSAB1", "KIT", "CPA3", "MS4A2", "TPSB2", "HPGDS"],
    "mural": ["ACTA2", "RGS5", "TAGLN", "MYH11", "MCAM", "DES"],
    "proliferating": ["MKI67", "TOP2A", "TYMS", "STMN1", "UBE2C", "PCLAF"],
}


class Tee(io.TextIOBase):
    def __init__(self, *streams: io.TextIOBase) -> None:
        self.streams = streams

    def write(self, text: str) -> int:
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


_log_handle = LOG_PATH.open("w", encoding="utf-8")
sys.stdout = Tee(sys.__stdout__, _log_handle)
sys.stderr = Tee(sys.__stderr__, _log_handle)


def log(msg: str) -> None:
    print(msg)


def normalize_gene_symbol(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().upper()
    text = re.sub(r"\s+", "", text)
    text = re.sub(r"///.*$", "", text)
    text = re.sub(r"//.*$", "", text)
    text = re.sub(r";.*$", "", text)
    text = re.sub(r",.*$", "", text)
    if text in {"", "NA", "N/A", "NULL", "---"}:
        return None
    return text


def read_frozen_modules(path: Path) -> dict[str, list[str]]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    modules = {}
    for entry in obj["modules"]:
        modules[entry["name"]] = [g for g in (normalize_gene_symbol(x) for x in entry["genes"]) if g]
    modules[obj["negative_control"]["name"]] = [g for g in (normalize_gene_symbol(x) for x in obj["negative_control"]["genes"]) if g]
    comp = obj.get("composite_definition")
    if comp:
        genes: list[str] = []
        for name in comp.get("includes", []):
            genes.extend(modules.get(name, []))
        modules[comp["name"]] = list(dict.fromkeys(genes))
    return modules


def read_csv_auto(path: Path, **kwargs: Any) -> pd.DataFrame:
    if path.suffix == ".gz":
        return pd.read_csv(path, compression="gzip", **kwargs)
    return pd.read_csv(path, **kwargs)


def collapse_duplicate_genes(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = [normalize_gene_symbol(x) for x in expr.index]
    expr = expr.loc[[x is not None for x in expr.index]]
    expr.index = pd.Index(expr.index)
    expr = expr.groupby(level=0).mean()
    return expr


def zscore_by_gene(expr: pd.DataFrame) -> pd.DataFrame:
    means = expr.mean(axis=1)
    sds = expr.std(axis=1, ddof=0).replace(0, np.nan)
    return expr.sub(means, axis=0).div(sds, axis=0)


def module_score(zmat: pd.DataFrame, genes: list[str]) -> pd.Series:
    genes_present = [g for g in genes if g in zmat.index]
    if not genes_present:
        return pd.Series(np.nan, index=zmat.columns)
    return zmat.loc[genes_present].mean(axis=0)


def standardize_series(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    sd = values.std(ddof=0)
    if pd.isna(sd) or sd == 0:
        return pd.Series(0.0, index=values.index)
    return (values - values.mean()) / sd


def hedges_g(case: pd.Series, ctrl: pd.Series) -> float:
    case = case.dropna().astype(float)
    ctrl = ctrl.dropna().astype(float)
    n1, n2 = len(case), len(ctrl)
    if n1 < 2 or n2 < 2:
        return np.nan
    sx = case.std(ddof=1)
    sy = ctrl.std(ddof=1)
    pooled = np.sqrt(((n1 - 1) * sx * sx + (n2 - 1) * sy * sy) / (n1 + n2 - 2))
    if pooled == 0 or pd.isna(pooled):
        return np.nan
    d = (case.mean() - ctrl.mean()) / pooled
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    return float(d * j)


def hedges_g_variance(g: float, n1: int, n2: int) -> float:
    if pd.isna(g) or n1 < 2 or n2 < 2:
        return np.nan
    return float((n1 + n2) / (n1 * n2) + (g * g) / (2 * (n1 + n2 - 2)))


def bootstrap_ci_stat(
    case: pd.Series,
    ctrl: pd.Series,
    n_boot: int = 4000,
    seed: int = 250704,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    x = case.dropna().astype(float).to_numpy()
    y = ctrl.dropna().astype(float).to_numpy()
    if len(x) < 2 or len(y) < 2:
        return (np.nan, np.nan)
    out = []
    for _ in range(n_boot):
        sx = rng.choice(x, size=len(x), replace=True)
        sy = rng.choice(y, size=len(y), replace=True)
        out.append(hedges_g(pd.Series(sx), pd.Series(sy)))
    return (float(np.nanpercentile(out, 2.5)), float(np.nanpercentile(out, 97.5)))


def exact_group_permutation(scores: pd.Series, labels: pd.Series) -> float:
    scores = pd.Series(scores, dtype=float).reset_index(drop=True)
    labels = pd.Series(labels, dtype=float).reset_index(drop=True)
    valid = scores.notna() & labels.notna()
    scores = scores.loc[valid].reset_index(drop=True)
    labels = labels.loc[valid].astype(int).reset_index(drop=True)
    n_case = int(labels.sum())
    n = len(scores)
    if n_case == 0 or n_case == n:
        return np.nan
    observed = abs(float(scores[labels == 1].mean() - scores[labels == 0].mean()))
    idx = list(range(n))
    extreme = 0
    total = 0
    for combo in itertools.combinations(idx, n_case):
        mask = np.zeros(n, dtype=bool)
        mask[list(combo)] = True
        diff = abs(float(scores[mask].mean() - scores[~mask].mean()))
        extreme += diff >= observed - 1e-12
        total += 1
    return float(extreme / total)


def paired_dz(diffs: pd.Series) -> float:
    diffs = diffs.dropna().astype(float)
    if len(diffs) < 2:
        return np.nan
    sd = diffs.std(ddof=1)
    if sd == 0 or pd.isna(sd):
        return np.nan
    return float(diffs.mean() / sd)


def matched_random_sets(expr: pd.DataFrame, target_genes: list[str], n_sets: int, seed: int = 250704) -> list[list[str]]:
    rng = np.random.default_rng(seed)
    mean_expr = expr.mean(axis=1)
    det = expr.notna().mean(axis=1)
    bins_expr = pd.qcut(mean_expr.rank(method="first"), q=min(20, mean_expr.shape[0]), labels=False, duplicates="drop")
    bins_det = pd.qcut(det.rank(method="first"), q=min(10, det.shape[0]), labels=False, duplicates="drop")
    annot = pd.DataFrame({"gene": expr.index.astype(str), "expr_bin": bins_expr, "det_bin": bins_det})
    gene_lookup = annot.set_index("gene")[["expr_bin", "det_bin"]]
    pair_pool = {key: grp["gene"].tolist() for key, grp in annot.groupby(["expr_bin", "det_bin"], dropna=False)}
    expr_pool = {key: grp["gene"].tolist() for key, grp in annot.groupby(["expr_bin"], dropna=False)}
    all_genes = annot["gene"].tolist()
    targets = [g for g in target_genes if g in set(all_genes)]
    out: list[list[str]] = []
    for _ in range(n_sets):
        used: set[str] = set()
        genes: list[str] = []
        for gene in targets:
            row = gene_lookup.loc[gene]
            pool = [g for g in pair_pool.get((row["expr_bin"], row["det_bin"]), []) if g not in used and g not in targets]
            if not pool:
                pool = [g for g in expr_pool.get(row["expr_bin"], []) if g not in used and g not in targets]
            if not pool:
                pool = [g for g in all_genes if g not in used and g not in targets]
            pick = str(rng.choice(pool))
            used.add(pick)
            genes.append(pick)
        out.append(genes)
    return out


def random_empirical_p(expr: pd.DataFrame, zmat: pd.DataFrame, target_genes: list[str], case_ids: list[str], ctrl_ids: list[str], observed_stat: float, n_sets: int) -> tuple[float, float]:
    null_stats = []
    for genes in matched_random_sets(expr, target_genes, n_sets=n_sets):
        score = module_score(zmat, genes)
        null_stats.append(hedges_g(score.loc[case_ids], score.loc[ctrl_ids]))
    null = pd.Series(null_stats, dtype=float)
    return float((null.abs() >= abs(observed_stat) - 1e-12).mean()), float((null < observed_stat).mean() * 100.0)


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pvals.astype(float)
    order = np.argsort(p.to_numpy())
    ranked = np.empty(len(p))
    vals = p.to_numpy()[order]
    m = len(vals)
    adj = np.minimum.accumulate((vals[::-1] * m / np.arange(m, 0, -1)))[::-1]
    ranked[order] = np.clip(adj, 0, 1)
    return pd.Series(ranked, index=p.index)


def random_effects_meta(df: pd.DataFrame, effect_col: str = "hedges_g", var_col: str = "variance") -> dict[str, Any]:
    work = df[[effect_col, var_col]].dropna().copy()
    if len(work) == 0:
        return {"k": 0, "pooled_effect": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "tau2": np.nan, "Q": np.nan, "I2": np.nan}
    yi = work[effect_col].to_numpy(dtype=float)
    vi = work[var_col].to_numpy(dtype=float)
    wi = 1.0 / vi
    mu_fixed = np.sum(wi * yi) / np.sum(wi)
    q = float(np.sum(wi * (yi - mu_fixed) ** 2))
    c = float(np.sum(wi) - (np.sum(wi ** 2) / np.sum(wi)))
    tau2 = max(0.0, (q - (len(yi) - 1)) / c) if c > 0 else 0.0
    wi_star = 1.0 / (vi + tau2)
    mu = float(np.sum(wi_star * yi) / np.sum(wi_star))
    se = float(np.sqrt(1.0 / np.sum(wi_star)))
    ci_low = mu - 1.96 * se
    ci_high = mu + 1.96 * se
    z = mu / se if se > 0 else np.nan
    pval = float(2 * stats.norm.sf(abs(z))) if pd.notna(z) else np.nan
    i2 = max(0.0, (q - (len(yi) - 1)) / q) * 100.0 if q > 0 else 0.0
    return {"k": int(len(yi)), "pooled_effect": mu, "ci_low": ci_low, "ci_high": ci_high, "se": se, "z": z, "p_value": pval, "tau2": tau2, "Q": q, "I2": i2}


def parse_probe_annotation(path: Path) -> pd.DataFrame:
    tbl = pd.read_csv(path, low_memory=False)
    probe_col = next((c for c in ["ID", "ID_REF", "Probe.Set.ID", "PROBEID"] if c in tbl.columns), None)
    symbol_col = next((c for c in ["Gene Symbol", "GENE_SYMBOL", "Gene symbol", "Symbol"] if c in tbl.columns), None)
    if probe_col is None or symbol_col is None:
        raise ValueError(f"Missing probe/gene annotation columns in {path}")
    ann = pd.DataFrame({"feature_id": tbl[probe_col].astype(str), "gene_symbol": tbl[symbol_col].map(normalize_gene_symbol)})
    ann = ann.dropna(subset=["gene_symbol"]).drop_duplicates("feature_id")
    return ann


def map_probes_to_genes(expr: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy().reset_index()
    expr = expr.rename(columns={expr.columns[0]: "feature_id"})
    expr["feature_id"] = expr["feature_id"].astype(str)
    merged = ann.merge(expr, on="feature_id", how="inner")
    if merged.empty:
        raise ValueError("No mapped probes retained")
    gene_expr = merged.drop(columns=["feature_id"]).set_index("gene_symbol")
    gene_expr = gene_expr.apply(pd.to_numeric, errors="coerce")
    gene_expr = gene_expr.groupby(level=0).mean()
    return gene_expr


def parse_group_from_text(text: str) -> str | None:
    t = str(text).lower()
    if "rheumatoid" in t or " ra " in f" {t} ":
        return "RA"
    if "osteoarthritis" in t or "osteoarthrisits" in t or "oateoarthritis" in t or "late stage oa" in t or "disease state: oa" in t:
        return "OA"
    if "healthy" in t or "normal" in t or "non-oa" in t or "control" in t:
        return "control"
    return None


def prepare_prior_processed_dataset(dataset_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(PRIOR_BULK / dataset_id / "expression.csv")
    expr = expr.set_index("gene_symbol")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = collapse_duplicate_genes(np.log2(expr + 1.0))
    meta = pd.read_csv(PRIOR_BULK / dataset_id / "metadata.csv")
    meta = meta.loc[meta["disease_group"].isin(["OA", "control"])].copy()
    expr = expr[meta["sample_id"].tolist()]
    meta["dataset_id"] = dataset_id
    meta["source_type"] = "prior_processed_gene_level"
    meta["control_context"] = meta["source_name"].fillna("").map(lambda x: "arthroscopy_or_possible_inflamed" if "biopsy" in str(x).lower() else "noninflamed_or_unspecified")
    return expr, meta


def prepare_gse12021() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(RAWGEO / "GSE12021" / "sample_metadata.csv", low_memory=False)
    meta["disease_group"] = meta.apply(
        lambda r: parse_group_from_text(" | ".join([str(r.get("title", "")), str(r.get("source_name_ch1", "")), str(r.get("characteristics_collapsed", ""))])),
        axis=1,
    )
    meta = meta.loc[meta["disease_group"].isin(["OA", "control", "RA"])].copy()
    parts = []
    for idx, gpl in [(1, "GPL96"), (2, "GPL97")]:
        expr = pd.read_csv(RAWGEO / "GSE12021" / f"series_{idx:02d}_expression_matrix.csv", low_memory=False).set_index("feature_id")
        expr = expr.apply(pd.to_numeric, errors="coerce")
        ann = parse_probe_annotation(PLATFORM_MAP[gpl])
        mapped = map_probes_to_genes(expr, ann)
        mapped = collapse_duplicate_genes(np.log2(mapped + 1.0))
        sub_meta = meta.loc[meta["platform"] == gpl].copy()
        common = [c for c in sub_meta["sample_id"].tolist() if c in mapped.columns]
        sub_meta = sub_meta.loc[sub_meta["sample_id"].isin(common)].copy()
        sub_meta["platform_block"] = gpl
        parts.append((mapped[common], sub_meta))
    expr_all = pd.concat([p[0] for p in parts], axis=1)
    meta_all = pd.concat([p[1] for p in parts], ignore_index=True)
    meta_all["dataset_id"] = "GSE12021"
    meta_all["source_type"] = "series_matrix_plus_platform_annotation"
    meta_all["control_context"] = "noninflamed_or_unspecified"
    return expr_all, meta_all


def prepare_gse82107() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(RAWGEO / "GSE82107" / "sample_metadata.csv", low_memory=False)
    meta["disease_group"] = meta.apply(
        lambda r: parse_group_from_text(" | ".join([str(r.get("title", "")), str(r.get("source_name_ch1", "")), str(r.get("characteristics_collapsed", "")), str(r.get("disease state:ch1", ""))])),
        axis=1,
    )
    meta = meta.loc[meta["disease_group"].isin(["OA", "control", "RA"])].copy()
    expr = pd.read_csv(RAWGEO / "GSE82107" / "series_01_expression_matrix.csv", low_memory=False).set_index("feature_id")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    ann = parse_probe_annotation(PLATFORM_MAP["GPL570"])
    mapped = map_probes_to_genes(expr, ann)
    mapped = collapse_duplicate_genes(np.log2(mapped + 1.0))
    common = [c for c in meta["sample_id"].tolist() if c in mapped.columns]
    meta = meta.loc[meta["sample_id"].isin(common)].copy()
    meta["dataset_id"] = "GSE82107"
    meta["source_type"] = "series_matrix_plus_platform_annotation"
    meta["control_context"] = "noninflamed_or_unspecified"
    return mapped[common], meta


def parse_dchip_xlsx(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name=0, header=None)
    header_row = None
    for i in range(min(50, len(frame))):
        if str(frame.iat[i, 0]).strip().lower() == "probe set":
            header_row = i
            break
    if header_row is None:
        raise ValueError("Could not locate dChip header row")
    header = frame.iloc[header_row].tolist()
    data = frame.iloc[header_row + 1 :].reset_index(drop=True)
    data.columns = header
    signal_cols = [c for c in data.columns[1:] if not str(c).endswith(" call")]
    expr = data[["probe set", *signal_cols]].copy()
    expr = expr.rename(columns={"probe set": "feature_id"}).set_index("feature_id")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.columns = [re.sub(r"\.cel$", "", str(c), flags=re.I) for c in expr.columns]
    return expr


def prepare_gse206848() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(RAWGEO / "GSE206848" / "sample_metadata.csv", low_memory=False)
    meta["disease_group"] = meta.apply(
        lambda r: parse_group_from_text(" | ".join([str(r.get("title", "")), str(r.get("source_name_ch1", "")), str(r.get("characteristics_collapsed", "")), str(r.get("tissue:ch1", ""))])),
        axis=1,
    )
    meta = meta.loc[meta["disease_group"].isin(["OA", "control", "RA"])].copy()
    meta["array_id"] = meta["supplementary_file"].astype(str).map(
        lambda x: re.sub(r"\.cel\.gz$", "", Path(x).name, flags=re.I).split("_", 1)[-1]
    )
    expr = parse_dchip_xlsx(RAWGEO / "GSE206848" / "GSE206848_synovium_All_expression.xlsx")
    ann = parse_probe_annotation(PLATFORM_MAP["GPL570"])
    mapped = map_probes_to_genes(expr, ann)
    mapped = collapse_duplicate_genes(np.log2(mapped + 1.0))
    common = [c for c in meta["array_id"].tolist() if c in mapped.columns]
    meta = meta.loc[meta["array_id"].isin(common)].copy()
    meta["sample_id"] = meta["array_id"]
    meta["dataset_id"] = "GSE206848"
    meta["source_type"] = "supplementary_dchip_expression"
    meta["control_context"] = meta["source_name_ch1"].fillna("").map(lambda x: "autopsy_nonoa" if "autopsy" in str(x).lower() or "normal synovium" in str(x).lower() else "noninflamed_or_unspecified")
    return mapped[common], meta


def prepare_gse89408() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(RAWGEO / "GSE89408" / "GSE89408_GEO_count_matrix_rename.txt.gz", sep="\t", compression="gzip", low_memory=False)
    first_col = expr.columns[0]
    expr = expr.rename(columns={first_col: "gene_symbol"}).set_index("gene_symbol")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    keep_cols = [c for c in expr.columns if c.startswith("normal_tissue_") or c.startswith("OA_tissue_")]
    expr = collapse_duplicate_genes(np.log2(expr[keep_cols] + 1.0))
    meta = pd.DataFrame({"sample_id": keep_cols})
    meta["disease_group"] = np.where(meta["sample_id"].str.startswith("OA_tissue_"), "OA", "control")
    meta["dataset_id"] = "GSE89408"
    meta["platform"] = "GPL11154"
    meta["source_type"] = "supplementary_count_matrix"
    meta["label_source"] = "supplementary matrix column names and GEO disease metadata"
    meta["control_context"] = "possibly_nonadult_or_heterogeneous_normal_biopsy"
    return expr[meta["sample_id"].tolist()], meta


def prepare_gse143514() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(RAWGEO / "GSE143514" / "GSE143514_mRNA_raw_count.txt.gz", sep="\t", compression="gzip", low_memory=False)
    expr = expr.rename(columns={expr.columns[0]: "gene_symbol"}).set_index("gene_symbol")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = collapse_duplicate_genes(np.log2(expr + 1.0))
    meta = pd.read_csv(RAWGEO / "GSE143514" / "sample_metadata.csv", low_memory=False)
    meta = meta.loc[meta["platform"] == "GPL20795"].copy()
    meta["matrix_name"] = meta["description"].astype(str).str.replace("_raw_count.txt", "", regex=False)
    meta["disease_group"] = meta["disease state:ch1"].map(lambda x: "OA" if str(x).strip().upper() == "OA" else ("control" if str(x).strip().lower() == "control" else None))
    meta = meta.loc[meta["disease_group"].isin(["OA", "control"])].copy()
    col_map = {
        "OA-1B": "Osteoarthritis_1B_count",
        "OA-2B": "Osteoarthritis_2B_count",
        "OA-3B": "Osteoarthritis_3B_count",
        "OA-5B": "Osteoarthritis_5B_count",
        "OA-8A": "Osteoarthritis_8A_count",
        "Nor-1B": "Normal_1B_count",
        "Nor-3B": "Normal_3B_count",
        "Nor-4A": "Normal_4A_count",
    }
    meta["sample_id"] = meta["matrix_name"].map(col_map)
    meta = meta.dropna(subset=["sample_id"]).copy()
    common = [c for c in meta["sample_id"].tolist() if c in expr.columns]
    meta = meta.loc[meta["sample_id"].isin(common)].copy()
    meta["dataset_id"] = "GSE143514"
    meta["source_type"] = "supplementary_mrna_raw_count"
    meta["control_context"] = "noninflamed_or_unspecified"
    return expr[common], meta


def prepare_legacy_dataset(dataset_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if dataset_id in {"GSE1919", "GSE55235", "GSE55457"}:
        return prepare_prior_processed_dataset(dataset_id)
    if dataset_id == "GSE12021":
        return prepare_gse12021()
    if dataset_id == "GSE82107":
        return prepare_gse82107()
    if dataset_id == "GSE206848":
        return prepare_gse206848()
    if dataset_id == "GSE89408":
        return prepare_gse89408()
    if dataset_id == "GSE143514":
        return prepare_gse143514()
    raise ValueError(dataset_id)


def dataset_module_effects(dataset_id: str, expr: pd.DataFrame, meta: pd.DataFrame, modules: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = meta.copy()
    if "platform_block" not in meta.columns:
        meta["platform_block"] = meta.get("platform", "single_platform")
    rows = []
    sample_audit_rows = []
    module_list = [PRIMARY_MODULE, *SECONDARY_MODULES, *SENSITIVITY_MODULES]
    for _, row in meta.iterrows():
        sample_audit_rows.append(
            {
                "dataset_id": dataset_id,
                "sample_id": row["sample_id"],
                "disease_group": row["disease_group"],
                "platform_block": row.get("platform_block", row.get("platform", "")),
                "label_source": row.get("label_source", "GEO sample metadata"),
                "control_context": row.get("control_context", "unspecified"),
                "source_type": row.get("source_type", ""),
            }
        )
    score_table = pd.DataFrame(index=meta["sample_id"].tolist())
    for block, sub_meta in meta.groupby("platform_block", dropna=False):
        sub_ids = sub_meta["sample_id"].tolist()
        zmat = zscore_by_gene(expr[sub_ids])
        for mod in module_list:
            sc = module_score(zmat, modules[mod]).loc[sub_ids]
            score_table.loc[sub_ids, mod] = standardize_series(sc).values
    score_table["disease_group"] = meta.set_index("sample_id")["disease_group"]
    for mod in module_list:
        case = score_table.loc[score_table["disease_group"] == "OA", mod]
        ctrl = score_table.loc[score_table["disease_group"] == "control", mod]
        g = hedges_g(case, ctrl)
        var = hedges_g_variance(g, len(case.dropna()), len(ctrl.dropna()))
        ci_low, ci_high = bootstrap_ci_stat(case, ctrl)
        rows.append(
            {
                "dataset_id": dataset_id,
                "module": mod,
                "n_oa": int(case.dropna().shape[0]),
                "n_control": int(ctrl.dropna().shape[0]),
                "hedges_g": g,
                "variance": var,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "direction": "OA_higher" if pd.notna(g) and g > 0 else ("OA_lower" if pd.notna(g) and g < 0 else "non_estimable"),
                "n_detected_genes": sum(gene in expr.index for gene in modules[mod]),
                "n_module_genes": len(modules[mod]),
                "pct_detected": 100.0 * sum(gene in expr.index for gene in modules[mod]) / len(modules[mod]),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(sample_audit_rows)


def run_legacy_meta(modules: dict[str, list[str]]) -> pd.DataFrame:
    log("Running legacy OA-versus-non-OA meta-analysis")
    datasets = ["GSE1919", "GSE12021", "GSE55235", "GSE55457", "GSE82107", "GSE206848", "GSE89408", "GSE143514"]
    qc_rows = []
    effect_rows = []
    audit_rows = []
    for dataset_id in datasets:
        expr, meta = prepare_legacy_dataset(dataset_id)
        qc_rows.append(
            {
                "dataset_id": dataset_id,
                "n_samples_total": int(meta.shape[0]),
                "n_oa": int((meta["disease_group"] == "OA").sum()),
                "n_control": int((meta["disease_group"] == "control").sum()),
                "n_ra": int((meta["disease_group"] == "RA").sum()) if "RA" in meta["disease_group"].values else 0,
                "n_genes": int(expr.shape[0]),
                "analysis_status": "included" if (meta["disease_group"] == "OA").sum() >= 2 and (meta["disease_group"] == "control").sum() >= 2 else "excluded",
                "source_type": meta["source_type"].iloc[0],
                "platforms": ";".join(sorted(meta["platform"].astype(str).unique())) if "platform" in meta.columns else "",
            }
        )
        eff, aud = dataset_module_effects(dataset_id, expr, meta.loc[meta["disease_group"].isin(["OA", "control"])].copy(), modules)
        effect_rows.append(eff)
        audit_rows.append(aud)
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(META1 / "dataset_qc.csv", index=False)
    sample_audit = pd.concat(audit_rows, ignore_index=True)
    sample_audit.to_csv(META1 / "sample_label_audit.csv", index=False)
    effects = pd.concat(effect_rows, ignore_index=True)
    effects.to_csv(META1 / "study_level_module_effects.csv", index=False)

    meta_rows = []
    for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
        sub = effects.loc[effects["module"] == mod].copy()
        re = random_effects_meta(sub)
        re["module"] = mod
        re["oa_positive_studies"] = int((sub["hedges_g"] > 0).sum())
        re["n_studies"] = int(sub["hedges_g"].notna().sum())
        meta_rows.append(re)
    meta_df = pd.DataFrame(meta_rows)
    secondary_mask = meta_df["module"].isin(SECONDARY_MODULES)
    meta_df.loc[secondary_mask, "bh_fdr"] = bh_fdr(meta_df.loc[secondary_mask, "p_value"])
    meta_df.to_csv(META1 / "random_effects_meta_results.csv", index=False)

    loo_rows = []
    mhcii = effects.loc[effects["module"] == PRIMARY_MODULE].copy()
    for ds in mhcii["dataset_id"].unique():
        re = random_effects_meta(mhcii.loc[mhcii["dataset_id"] != ds])
        loo_rows.append({"left_out_study": ds, **re})
    loo_df = pd.DataFrame(loo_rows)
    loo_df.to_csv(META1 / "leave_one_study_out.csv", index=False)

    hetero = meta_df[["module", "k", "Q", "I2", "tau2"]].copy()
    hetero.to_csv(META1 / "heterogeneity_results.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    plot_df = mhcii[["dataset_id", "hedges_g", "ci_low", "ci_high"]].copy()
    plot_df = plot_df.sort_values("hedges_g")
    y = np.arange(len(plot_df))
    ax.errorbar(plot_df["hedges_g"], y, xerr=[plot_df["hedges_g"] - plot_df["ci_low"], plot_df["ci_high"] - plot_df["hedges_g"]], fmt="o", color="#4c78a8", ecolor="#4c78a8", capsize=3)
    pooled = meta_df.loc[meta_df["module"] == PRIMARY_MODULE].iloc[0]
    ax.axvline(0.0, color="0.5", linestyle="--", linewidth=1)
    ax.axvspan(pooled["ci_low"], pooled["ci_high"], color="#f28e2b", alpha=0.2)
    ax.axvline(pooled["pooled_effect"], color="#f28e2b", linewidth=2)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["dataset_id"])
    ax.set_xlabel("Hedges' g (OA vs non-OA)")
    ax.set_title("Legacy synovium MHC-II meta-analysis")
    fig.tight_layout()
    fig.savefig(META1 / "legacy_meta_forest_plot.png", dpi=180)
    plt.close(fig)

    sensitivity_rows = []
    exclude_contexts = {"arthroscopy_or_possible_inflamed", "possibly_nonadult_or_heterogeneous_normal_biopsy"}
    eligible = sample_audit.loc[~sample_audit["control_context"].isin(exclude_contexts), "dataset_id"].unique().tolist()
    re = random_effects_meta(mhcii.loc[mhcii["dataset_id"].isin(eligible)])
    sensitivity_rows.append({"analysis": "exclude_arthroscopy_or_potentially_inflamed_controls", **re})
    pd.DataFrame(sensitivity_rows).to_csv(META1 / "heterogeneity_results_sensitivity_tmp.csv", index=False)

    interp = [
        "# Legacy Meta Interpretation",
        "",
        f"- Primary pooled MHC-II effect: {pooled['pooled_effect']:.3f} (95% CI {pooled['ci_low']:.3f} to {pooled['ci_high']:.3f}), p={pooled['p_value']:.4g}.",
        f"- OA-positive study direction count: {int(pooled['oa_positive_studies'])} of {int(pooled['n_studies'])}.",
        f"- Leave-one-study-out pooled direction remained {'positive' if (loo_df['pooled_effect'] > 0).all() else 'mixed'} across the legacy panel.",
        "- Secondary frozen modules were meta-analyzed separately and FDR-adjusted at the module level.",
    ]
    (META1 / "legacy_meta_interpretation.md").write_text("\n".join(interp) + "\n", encoding="utf-8")
    return meta_df


def parse_gse283079_metadata() -> pd.DataFrame:
    meta = pd.read_csv(RAWGEO / "GSE283079" / "sample_metadata.csv", low_memory=False)
    meta["sample_matrix_id"] = meta["title"].astype(str).str.replace("Knee ", "", regex=False).str.replace(" ", "_", regex=False)
    meta["sample_matrix_id"] = meta["sample_matrix_id"].str.replace("osteoarthritis_", "OA_", regex=False)
    meta["sample_matrix_id"] = meta["sample_matrix_id"].str.replace("non-osteoarthritis_", "non-OA_", regex=False)
    meta["sample_matrix_id"] = meta["sample_matrix_id"].str.replace("synovium_", "syn_", regex=False)
    meta["disease_group"] = meta["disease state:ch1"].map(lambda x: "OA" if "osteoarthritis" in str(x).lower() and "non-" not in str(x).lower() else ("control" if "non-osteoarthritis" in str(x).lower() else None))
    meta["control_context"] = meta["disease_group"].map(lambda x: "acl_reconstruction_arthroscopy_control" if x == "control" else "oa_tka")
    return meta.loc[meta["disease_group"].isin(["OA", "control"])].copy()


def run_gse283079_holdout(modules: dict[str, list[str]]) -> pd.DataFrame:
    log("Running GSE283079 hold-out analysis")
    meta = parse_gse283079_metadata()
    meta.to_csv(HOLD2 / "sample_label_audit.csv", index=False)
    expr = pd.read_csv(RAWGEO / "GSE283079" / "GSE283079_gene_count_RNAseq.txt.gz", sep="\t", compression="gzip", low_memory=False)
    expr = expr.rename(columns={expr.columns[0]: "gene_symbol"}).set_index("gene_symbol")
    expr.index = expr.index.map(lambda x: normalize_gene_symbol(str(x).split("_", 1)[1] if "_" in str(x) else x))
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[expr.mean(axis=1) > 1].copy()
    expr = collapse_duplicate_genes(np.log2(expr + 1.0))
    common = [c for c in meta["sample_matrix_id"].tolist() if c in expr.columns]
    meta = meta.loc[meta["sample_matrix_id"].isin(common)].copy()
    expr = expr[common]
    zmat = zscore_by_gene(expr)

    rows = []
    for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
        score = standardize_series(module_score(zmat, modules[mod]))
        case_ids = meta.loc[meta["disease_group"] == "OA", "sample_matrix_id"].tolist()
        ctrl_ids = meta.loc[meta["disease_group"] == "control", "sample_matrix_id"].tolist()
        g = hedges_g(score.loc[case_ids], score.loc[ctrl_ids])
        ci_low, ci_high = bootstrap_ci_stat(score.loc[case_ids], score.loc[ctrl_ids])
        p_exact = exact_group_permutation(score.loc[case_ids + ctrl_ids], pd.Series([1] * len(case_ids) + [0] * len(ctrl_ids)))
        emp_p, perc = random_empirical_p(expr, zmat, modules[mod], case_ids, ctrl_ids, g, n_sets=3000 if mod == PRIMARY_MODULE else 1000)
        rows.append({"module": mod, "n_oa": len(case_ids), "n_control": len(ctrl_ids), "hedges_g": g, "ci_low": ci_low, "ci_high": ci_high, "exact_permutation_p": p_exact, "random_gene_empirical_p": emp_p, "random_gene_percentile": perc})
    results = pd.DataFrame(rows)
    results.to_csv(HOLD2 / "module_score_results.csv", index=False)

    mhcii = results.loc[results["module"] == PRIMARY_MODULE].iloc[0]
    generic_score = standardize_series(module_score(zmat, GENERIC_MHCII_REFERENCE))
    hk_score = standardize_series(module_score(zmat, HOUSEKEEPING_CONTROL))
    case_ids = meta.loc[meta["disease_group"] == "OA", "sample_matrix_id"].tolist()
    ctrl_ids = meta.loc[meta["disease_group"] == "control", "sample_matrix_id"].tolist()
    generic_g = hedges_g(generic_score.loc[case_ids], generic_score.loc[ctrl_ids])
    hk_g = hedges_g(hk_score.loc[case_ids], hk_score.loc[ctrl_ids])
    random_df = pd.DataFrame(
        [
            {"control_type": "expression_matched_random_gene_sets", "effect_size": mhcii["hedges_g"], "empirical_p": mhcii["random_gene_empirical_p"]},
            {"control_type": "generic_mhcii_reference", "effect_size": generic_g},
            {"control_type": "housekeeping_like_control", "effect_size": hk_g},
        ]
    )
    random_df.to_csv(HOLD2 / "random_gene_set_results.csv", index=False)

    logo_rows = []
    target_genes = modules[PRIMARY_MODULE]
    for gene in target_genes:
        score = standardize_series(module_score(zmat, [g for g in target_genes if g != gene]))
        logo_rows.append({"left_out_gene": gene, "hedges_g": hedges_g(score.loc[case_ids], score.loc[ctrl_ids])})
    pd.DataFrame(logo_rows).to_csv(HOLD2 / "leave_one_gene_out.csv", index=False)

    text = [
        "# GSE283079 Hold-out Interpretation",
        "",
        "- This dataset serves as a fixed external OA-versus-non-OA disease-state hold-out.",
        "- Non-OA controls were obtained during arthroscopy after ACL reconstruction and are not healthy volunteer controls.",
        f"- Primary MHC-II effect: {mhcii['hedges_g']:.3f} (95% CI {mhcii['ci_low']:.3f} to {mhcii['ci_high']:.3f}), exact permutation p={mhcii['exact_permutation_p']:.4g}.",
        f"- Generic MHC-II reference effect: {generic_g:.3f}; housekeeping-like control effect: {hk_g:.3f}.",
    ]
    (HOLD2 / "holdout_interpretation.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    return results


def read_tar_member_text(tf: tarfile.TarFile, member_name: str) -> list[str]:
    fh = tf.extractfile(member_name)
    assert fh is not None
    with gzip.open(fh, "rt", encoding="utf-8", errors="ignore") as handle:
        return [line.rstrip("\n") for line in handle]


def read_scrna_sample(tf: tarfile.TarFile, sample_prefix: str) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    members = [m.name for m in tf.getmembers()]
    mtx_member = next(m for m in members if m.endswith(f"{sample_prefix}_matrix.mtx.gz"))
    feat_member = next(m for m in members if m.endswith(f"{sample_prefix}_features.tsv.gz"))
    bar_member = next(m for m in members if m.endswith(f"{sample_prefix}_barcodes.tsv.gz"))
    fh = tf.extractfile(mtx_member)
    assert fh is not None
    mat = mmread(gzip.GzipFile(fileobj=fh)).tocsr()
    genes = [normalize_gene_symbol(line.split("\t")[1] if "\t" in line else line) for line in read_tar_member_text(tf, feat_member)]
    cells = [line.strip() for line in read_tar_member_text(tf, bar_member)]
    keep = [i for i, g in enumerate(genes) if g is not None]
    mat = mat[keep, :]
    genes = [genes[i] for i in keep]
    return mat, genes, cells


def sparse_log1p_cpm(mat: sparse.csr_matrix) -> sparse.csr_matrix:
    mat = mat.tocsc(copy=True).astype(np.float32)
    lib = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)
    lib[lib == 0] = 1.0
    scale = 1e6 / lib
    repeats = np.diff(mat.indptr)
    mat.data *= np.repeat(scale, repeats)
    mat.data = np.log1p(mat.data)
    return mat.tocsr()


def score_marker_sets(norm: sparse.csr_matrix, genes: list[str]) -> pd.DataFrame:
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    rows = []
    for celltype, markers in CELLTYPE_MARKERS.items():
        idx = [gene_to_idx[g] for g in markers if g in gene_to_idx]
        if not idx:
            rows.append(pd.Series(np.nan, index=np.arange(norm.shape[1]), name=celltype))
        else:
            sub = norm[idx, :]
            rows.append(pd.Series(np.asarray(sub.mean(axis=0)).ravel(), name=celltype))
    return pd.DataFrame(rows)


def annotate_cells(marker_scores: pd.DataFrame) -> pd.DataFrame:
    assign = marker_scores.idxmax(axis=0)
    ordered = np.sort(marker_scores.to_numpy(), axis=0)
    top = ordered[-1, :]
    second = ordered[-2, :] if marker_scores.shape[0] > 1 else np.repeat(np.nan, marker_scores.shape[1])
    return pd.DataFrame({"assigned_cell_type": assign.values, "top_score": top, "second_score": second})


def pseudobulk_counts(mat: sparse.csr_matrix, labels: pd.Series, genes: list[str], donor_id: str) -> pd.DataFrame:
    rows = []
    cols = []
    for label in sorted(labels.dropna().unique().tolist()):
        idx = np.where(labels.values == label)[0]
        if len(idx) == 0:
            continue
        vec = np.asarray(mat[:, idx].sum(axis=1)).ravel()
        rows.append(vec)
        cols.append(f"{donor_id}|{label}")
    if not rows:
        return pd.DataFrame(index=genes)
    out = pd.DataFrame(np.vstack(rows).T, index=genes, columns=cols)
    return out.groupby(level=0).sum()


def run_gse216651_localization(modules: dict[str, list[str]]) -> pd.DataFrame:
    log("Running GSE216651 localization analysis")
    meta = pd.read_csv(RAWGEO / "GSE216651" / "sample_metadata.csv", low_memory=False)
    meta["donor_id"] = meta["sample_id"]
    meta["disease_group"] = meta["disease state:ch1"].map(lambda x: "OA" if "osteoarthritis" in str(x).lower() else ("healthy" if "healthy" in str(x).lower() else None))
    meta["sample_modality"] = meta["title"].astype(str).map(lambda x: "sc" if x.startswith("sc_") else ("sn" if x.startswith("sn_") else "unknown"))
    meta["tissue_restriction_status"] = "joint_ipfp_synovium_combined_not_separable"
    meta[["donor_id", "sample_id", "title", "disease_group", "sample_modality", "tissue_restriction_status"]].to_csv(LOC3 / "donor_metadata.csv", index=False)

    ann_rows = []
    pb_list = []
    with tarfile.open(RAWGEO / "GSE216651" / "GSE216651_RAW.tar") as tf:
        for _, row in meta.iterrows():
            sample_prefix = row["title"]
            mat, genes, cells = read_scrna_sample(tf, sample_prefix)
            norm = sparse_log1p_cpm(mat)
            marker_scores = score_marker_sets(norm, genes)
            annot = annotate_cells(marker_scores)
            annot["cell_id"] = cells
            annot["donor_id"] = row["donor_id"]
            ann_rows.append(annot[["cell_id", "donor_id", "assigned_cell_type", "top_score", "second_score"]])
            pb_list.append(pseudobulk_counts(mat, annot["assigned_cell_type"], genes, row["donor_id"]))
    ann_df = pd.concat(ann_rows, ignore_index=True)
    ann_df.to_csv(LOC3 / "cell_annotation_audit.csv", index=False)

    pb_counts = pd.concat(pb_list, axis=1).fillna(0.0)
    pb_counts = pb_counts.groupby(level=0).sum()
    pb_norm = np.log1p(pb_counts.div(pb_counts.sum(axis=0), axis=1) * 1e6)
    pb_norm = collapse_duplicate_genes(pb_norm)
    pb_zmat = zscore_by_gene(pb_norm)

    score_rows = []
    for col in pb_norm.columns:
        donor, celltype = col.split("|", 1)
        row = {"donor_id": donor, "cell_type": celltype, "n_cells": int((ann_df["donor_id"].eq(donor) & ann_df["assigned_cell_type"].eq(celltype)).sum())}
        for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
            row[mod] = float(module_score(pb_zmat[[col]], modules[mod]).iloc[0])
        score_rows.append(row)
    pb_scores = pd.DataFrame(score_rows).merge(meta[["donor_id", "disease_group", "sample_modality", "tissue_restriction_status"]], on="donor_id", how="left")
    pb_scores.to_csv(LOC3 / "donor_pseudobulk_scores.csv", index=False)

    loc_rows = []
    for donor in sorted(pb_scores["donor_id"].unique()):
        sub = pb_scores.loc[pb_scores["donor_id"] == donor].copy()
        if "myeloid" not in sub["cell_type"].tolist():
            continue
        nonmy = sub.loc[sub["cell_type"] != "myeloid"]
        if nonmy.empty:
            continue
        for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
            loc_rows.append(
                {
                    "donor_id": donor,
                    "module": mod,
                    "myeloid_score": sub.loc[sub["cell_type"] == "myeloid", mod].iloc[0],
                    "nonmyeloid_mean_score": float(nonmy[mod].mean()),
                    "difference": float(sub.loc[sub["cell_type"] == "myeloid", mod].iloc[0] - nonmy[mod].mean()),
                }
            )
    loc_df = pd.DataFrame(loc_rows)
    myeloid = pb_scores.loc[pb_scores["cell_type"] == "myeloid"].copy()
    case_ids = myeloid.loc[myeloid["disease_group"] == "OA", "donor_id"].tolist()
    ctrl_ids = myeloid.loc[myeloid["disease_group"] == "healthy", "donor_id"].tolist()
    expr_my = pb_norm[[f"{d}|myeloid" for d in myeloid["donor_id"]]]
    expr_my.columns = myeloid["donor_id"].tolist()
    zmat_my = zscore_by_gene(expr_my)
    out_rows = []
    for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
        score = myeloid.set_index("donor_id")[mod]
        g = hedges_g(score.loc[case_ids], score.loc[ctrl_ids])
        ci_low, ci_high = bootstrap_ci_stat(score.loc[case_ids], score.loc[ctrl_ids])
        emp_p, perc = random_empirical_p(expr_my, zmat_my, modules[mod], case_ids, ctrl_ids, g, n_sets=1000 if mod == PRIMARY_MODULE else 500)
        out_rows.append(
            {
                "module": mod,
                "oa_vs_healthy_hedges_g": g,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "expression_matched_random_empirical_p": emp_p,
                "expression_matched_random_percentile": perc,
                "mean_myeloid_minus_nonmyeloid": float(loc_df.loc[loc_df["module"] == mod, "difference"].mean()) if not loc_df.empty else np.nan,
                "n_oa_donors": len(case_ids),
                "n_healthy_donors": len(ctrl_ids),
                "synovium_only_primary_estimable": "no",
            }
        )
    loc_results = pd.DataFrame(out_rows)
    loc_results.to_csv(LOC3 / "myeloid_localization_results.csv", index=False)

    lodo_rows = []
    for donor in myeloid["donor_id"]:
        sub = myeloid.loc[myeloid["donor_id"] != donor].copy()
        case = sub.loc[sub["disease_group"] == "OA", PRIMARY_MODULE]
        ctrl = sub.loc[sub["disease_group"] == "healthy", PRIMARY_MODULE]
        lodo_rows.append({"left_out_donor": donor, "hedges_g": hedges_g(case, ctrl)})
    pd.DataFrame(lodo_rows).to_csv(LOC3 / "leave_one_donor_out.csv", index=False)

    mhcii = loc_results.loc[loc_results["module"] == PRIMARY_MODULE].iloc[0]
    generic_score = standardize_series(module_score(zmat_my, GENERIC_MHCII_REFERENCE))
    hk_score = standardize_series(module_score(zmat_my, HOUSEKEEPING_CONTROL))
    pd.DataFrame(
        [
            {"control_type": "expression_matched_random_gene_sets", "effect_size": mhcii["oa_vs_healthy_hedges_g"], "empirical_p": mhcii["expression_matched_random_empirical_p"]},
            {"control_type": "generic_mhcii_reference", "effect_size": hedges_g(generic_score.loc[case_ids], generic_score.loc[ctrl_ids])},
            {"control_type": "housekeeping_like_control", "effect_size": hedges_g(hk_score.loc[case_ids], hk_score.loc[ctrl_ids])},
        ]
    ).to_csv(LOC3 / "random_gene_set_results.csv", index=False)

    interp = [
        "# GSE216651 Localization Interpretation",
        "",
        "- Public processed files support donor-level OA-versus-healthy joint-cell pseudobulk analysis and myeloid localization assessment.",
        "- Synovium-only primary restriction was not estimable because the public donor-level labels combine IPFP and synovium rather than separating them.",
        "- Cell-level scores were used for annotation and visualization only; donor-level pseudobulk remained the unit of inference.",
    ]
    (LOC3 / "localization_interpretation.md").write_text("\n".join(interp) + "\n", encoding="utf-8")
    return loc_results


def run_gse283080_exploratory(modules: dict[str, list[str]]) -> pd.DataFrame:
    log("Running GSE283080 exploratory cellular-state localization")
    meta = pd.read_csv(RAWGEO / "GSE283080" / "sample_metadata.csv", low_memory=False)
    meta["state_group"] = meta["title"].astype(str).map(lambda x: "inflammatory" if "infla" in x.lower() else ("fibrotic" if "fibro" in x.lower() else None))
    meta["donor_id"] = meta["sample_id"]
    meta = meta.loc[meta["state_group"].isin(["inflammatory", "fibrotic"])].copy()
    meta[["donor_id", "sample_id", "title", "state_group"]].to_csv(EXP4 / "sample_metadata_audit.csv", index=False)

    ann_rows = []
    pb_list = []
    with tarfile.open(RAWGEO / "GSE283080" / "GSE283080_RAW.tar") as tf:
        members = [m.name for m in tf.getmembers()]
        mapping = {
            "GSM8655601": "Infla1",
            "GSM8655602": "Infla2",
            "GSM8655603": "Fibro1",
            "GSM8655604": "Fibro2",
        }
        for _, row in meta.iterrows():
            prefix = mapping[row["sample_id"]]
            mtx_member = next(m for m in members if m.endswith(f"{prefix}_matrix.mtx.gz"))
            feat_member = next(m for m in members if m.endswith(f"{prefix}_features.tsv.gz"))
            bar_member = next(m for m in members if m.endswith(f"{prefix}_barcodes.tsv.gz"))
            fh = tf.extractfile(mtx_member)
            assert fh is not None
            mat = mmread(gzip.GzipFile(fileobj=fh)).tocsr()
            genes = [normalize_gene_symbol(line.split("\t")[1] if "\t" in line else line) for line in read_tar_member_text(tf, feat_member)]
            cells = [line.strip() for line in read_tar_member_text(tf, bar_member)]
            keep = [i for i, g in enumerate(genes) if g is not None]
            mat = mat[keep, :]
            genes = [genes[i] for i in keep]
            norm = sparse_log1p_cpm(mat)
            marker_scores = score_marker_sets(norm, genes)
            annot = annotate_cells(marker_scores)
            annot["cell_id"] = cells
            annot["donor_id"] = row["donor_id"]
            ann_rows.append(annot)
            pb_list.append(pseudobulk_counts(mat, annot["assigned_cell_type"], genes, row["donor_id"]))

    ann_df = pd.concat(ann_rows, ignore_index=True)
    pb_counts = pd.concat(pb_list, axis=1).fillna(0.0)
    pb_counts = pb_counts.groupby(level=0).sum()
    pb_norm = np.log1p(pb_counts.div(pb_counts.sum(axis=0), axis=1) * 1e6)
    pb_norm = collapse_duplicate_genes(pb_norm)
    pb_zmat = zscore_by_gene(pb_norm)

    rows = []
    for col in pb_norm.columns:
        donor, celltype = col.split("|", 1)
        row = {"donor_id": donor, "cell_type": celltype, "state_group": meta.set_index("donor_id").loc[donor, "state_group"]}
        for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
            row[mod] = float(module_score(pb_zmat[[col]], modules[mod]).iloc[0])
        rows.append(row)
    cell_scores = pd.DataFrame(rows)
    cell_scores.to_csv(EXP4 / "cell_state_module_scores.csv", index=False)

    dir_rows = []
    for mod in [PRIMARY_MODULE, *SECONDARY_MODULES]:
        for celltype in sorted(cell_scores["cell_type"].unique()):
            sub = cell_scores.loc[cell_scores["cell_type"] == celltype].copy()
            case = sub.loc[sub["state_group"] == "inflammatory", mod]
            ctrl = sub.loc[sub["state_group"] == "fibrotic", mod]
            dir_rows.append(
                {
                    "module": mod,
                    "cell_type": celltype,
                    "direction": "inflammatory_higher" if case.mean() > ctrl.mean() else "fibrotic_higher",
                    "mean_difference": float(case.mean() - ctrl.mean()),
                    "n_inflammatory": int(case.shape[0]),
                    "n_fibrotic": int(ctrl.shape[0]),
                }
            )
    pd.DataFrame(dir_rows).to_csv(EXP4 / "directionality_summary.csv", index=False)
    text = [
        "# GSE283080 Exploratory Interpretation",
        "",
        "- This dataset is restricted to exploratory cellular-state localization.",
        "- Sample size is n=2 inflammatory versus n=2 fibroblast/fibrotic OA donors.",
        "- No confirmatory p values were used as evidence.",
    ]
    (EXP4 / "exploratory_interpretation.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    return cell_scores


def synthesize(meta_df: pd.DataFrame, holdout_df: pd.DataFrame, loc_df: pd.DataFrame) -> None:
    log("Synthesizing cross-layer evidence")
    symptom_996 = pd.read_csv(PRIOR_SYMPTOM / "02_gse99662_audit" / "reproducibility_results.csv")
    symptom_176 = pd.read_csv(PRIOR_SYMPTOM / "03_gse176223_audit" / "strict_paired_results.csv")
    symptom_248 = pd.read_csv(PRIOR_SYMPTOM / "04_independent_scrna" / "phenotype_association_results.csv")

    legacy_primary = meta_df.loc[meta_df["module"] == PRIMARY_MODULE].iloc[0]
    hold_primary = holdout_df.loc[holdout_df["module"] == PRIMARY_MODULE].iloc[0]
    loc_primary = loc_df.loc[loc_df["module"] == PRIMARY_MODULE].iloc[0]
    g996 = symptom_996.loc[symptom_996["module"] == PRIMARY_MODULE].iloc[0]
    g176 = symptom_176.loc[symptom_176["module"] == PRIMARY_MODULE].iloc[0]
    g248 = symptom_248.loc[symptom_248["module"] == PRIMARY_MODULE].iloc[0]

    final_summary = pd.DataFrame(
        [
            {"layer": "disease_state_conservation", "dataset_or_analysis": "legacy_meta_mhcii", "effect": legacy_primary["pooled_effect"], "ci_low": legacy_primary["ci_low"], "ci_high": legacy_primary["ci_high"], "interpretation_class": "confirmatory"},
            {"layer": "disease_state_conservation", "dataset_or_analysis": "gse283079_holdout_mhcii", "effect": hold_primary["hedges_g"], "ci_low": hold_primary["ci_low"], "ci_high": hold_primary["ci_high"], "interpretation_class": "holdout"},
            {"layer": "cellular_localization", "dataset_or_analysis": "gse216651_myeloid_mhcii", "effect": loc_primary["oa_vs_healthy_hedges_g"], "ci_low": loc_primary["ci_low"], "ci_high": loc_primary["ci_high"], "interpretation_class": "localization"},
            {"layer": "symptom_state", "dataset_or_analysis": "gse99662_high_vs_low_pain_mhcii", "effect": g996["hedges_g"], "ci_low": g996["ci_low"], "ci_high": g996["ci_high"], "interpretation_class": "non_supportive"},
            {"layer": "symptom_state", "dataset_or_analysis": "gse176223_paired_pain_vs_no_pain_mhcii", "effect": g176["paired_dz"], "ci_low": g176["ci_low"], "ci_high": g176["ci_high"], "interpretation_class": "exploratory"},
            {"layer": "symptom_state", "dataset_or_analysis": "gse248453_myeloid_high_vs_low_pain_mhcii", "effect": g248["hedges_g"], "ci_low": g248["ci_low"], "ci_high": g248["ci_high"], "interpretation_class": "non_supportive"},
        ]
    )
    final_summary.to_csv(SYN5 / "final_cross_layer_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    plot_df = final_summary.copy()
    y = np.arange(len(plot_df))
    ax.errorbar(plot_df["effect"], y, xerr=[plot_df["effect"] - plot_df["ci_low"], plot_df["ci_high"] - plot_df["effect"]], fmt="o", color="#4c78a8", ecolor="#4c78a8", capsize=3)
    ax.axvline(0.0, color="0.5", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["dataset_or_analysis"])
    ax.set_xlabel("Effect size")
    ax.set_title("Cross-layer effect summary")
    fig.tight_layout()
    fig.savefig(SYN5 / "final_effect_size_plot.png", dpi=180)
    plt.close(fig)

    legacy_supported = (legacy_primary["pooled_effect"] > 0) and (legacy_primary["ci_low"] > 0) and (legacy_primary["oa_positive_studies"] >= 6)
    holdout_supported = (hold_primary["hedges_g"] > 0)
    localization_supported = loc_primary["mean_myeloid_minus_nonmyeloid"] > 0 if pd.notna(loc_primary["mean_myeloid_minus_nonmyeloid"]) else False
    symptom_not_universal = (g996["hedges_g"] < 0) and (g248["hedges_g"] < 0)

    if legacy_supported and holdout_supported and localization_supported and symptom_not_universal:
        decision = "JTM-oriented Conditional Go"
        next_step = "independent tissue validation with pain and/or MRI synovitis metadata"
    elif legacy_supported:
        decision = "OAC-oriented Reframe"
        next_step = "OA-native local myeloid tissue-state manuscript or tissue validation"
    else:
        decision = "No-Go for JTM rescue"
        next_step = "terminate the JIA-derived MHC-II rescue attempt"

    report = [
        "# Final Decision Report",
        "",
        f"- Final decision class: **{decision}**",
        "",
        "## Confirmed",
        f"- Legacy OA-versus-non-OA disease-state MHC-II conservation is {'supported' if legacy_supported else 'not supported'} by the locked meta-analysis panel.",
        f"- Existing symptom-state evidence continues to show that MHC-II is not a universal OA pain marker.",
        "",
        "## Exploratory",
        "- S100A8/A9 and complement/coagulation remain exploratory local pain-site signals only.",
        "- GSE283080 is exploratory cellular-state visualization only.",
        "",
        "## Unsupported or bounded",
        "- No direct transcriptomic validation claim is made in OAI or FNIH.",
        "- GSE216651 did not allow strict synovium-only restriction because public processed labels combine IPFP and synovium.",
        "",
        "## Next permissible action",
        f"- {next_step}",
    ]
    (SYN5 / "final_decision_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> int:
    log("Loading frozen modules")
    modules = read_frozen_modules(FROZEN_YAML)
    legacy_meta = run_legacy_meta(modules)
    holdout = run_gse283079_holdout(modules)
    loc = run_gse216651_localization(modules)
    run_gse283080_exploratory(modules)
    synthesize(legacy_meta, holdout, loc)
    log("Layer 2.5S completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
