#!/usr/bin/env python3
"""Compile separated YPT statistics plex exports into website data assets."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd


SOURCE_DATASETS = {
    "os": ("os_0.7", "One-Shot (Complete)"),
    "frac": ("frac_0.7", "Fractionated"),
}
META_COLUMNS = ["Protein Id", "gene_symbol", "prot_description", "Site Position", "sequence"]
PART_SIZE = 8000
HIT_R = 2.0
SIG_NEGLOG10 = -math.log10(0.05)
PVALUE_FLOOR = 1e-300


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, separators=(",", ":"), allow_nan=False)


def clean_sequence(value) -> str:
    text = str(value or "")
    parts = text.split(".")
    if len(parts) >= 3:
        text = parts[1]
    return re.sub(r"[^A-Za-z]", "", text)


def clean_site(value) -> str:
    text = str(value or "").strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def site_positions(site: str) -> list[str]:
    found = re.findall(r"\d+", site)
    return found or [site]


def uniprot_from_protein_id(value: str) -> str:
    parts = str(value or "").split("|")
    return parts[1] if len(parts) >= 3 else str(value or "")


def site_label(gene: str, site: str) -> str:
    suffix = "_".join(site_positions(site))
    return f"{gene}_Y{suffix}" if suffix else gene


def metadata_keys(row: dict) -> list[tuple[str, ...]]:
    protein = str(row.get("proteinId") or "")
    gene = str(row.get("gene") or "")
    site = clean_site(row.get("site") or "")
    seq = clean_sequence(row.get("sequence") or "")
    label = str(row.get("label") or "")
    return [
        (protein, gene, site, seq),
        (gene, site, seq),
        (protein, gene, site),
        (label,),
    ]


def build_existing_metadata(asset_root: Path, manifest: dict) -> dict[tuple[str, ...], dict]:
    mapping: dict[tuple[str, ...], dict] = {}
    for dataset in manifest.get("datasets", {}).values():
        catalog_path = asset_root / dataset.get("catalog", "")
        if not catalog_path.exists():
            continue
        for row in read_json(catalog_path):
            meta = {
                "contactTypes": row.get("contactTypes", []),
                "contactKeys": row.get("contactKeys", []),
                "x": row.get("x"),
                "y": row.get("y"),
            }
            for key in metadata_keys(row):
                mapping.setdefault(key, meta)
    return mapping


def source_row_key(row) -> tuple[str, str, str, str]:
    protein = str(row["Protein Id"])
    gene = str(row["gene_symbol"])
    site = clean_site(row["Site Position"])
    seq = clean_sequence(row["sequence"])
    return protein, gene, site, seq


def source_metadata(row, existing_meta: dict[tuple[str, ...], dict]) -> dict:
    protein, gene, site, seq = source_row_key(row)
    label = site_label(gene, site)
    lookup_keys = [
        (protein, gene, site, seq),
        (gene, site, seq),
        (protein, gene, site),
        (label,),
    ]
    matched = {}
    for key in lookup_keys:
        if key in existing_meta:
            matched = existing_meta[key]
            break
    return {
        "label": label,
        "gene": gene,
        "site": site,
        "positions": site_positions(site),
        "proteinId": protein,
        "uniprot": uniprot_from_protein_id(protein),
        "description": str(row.get("prot_description") or ""),
        "sequence": seq,
        "contactTypes": matched.get("contactTypes", []),
        "contactKeys": matched.get("contactKeys", []),
        "x": matched.get("x"),
        "y": matched.get("y"),
    }


def finite_number(value, fallback=None):
    try:
        num = float(value)
    except (TypeError, ValueError):
        return fallback
    return num if math.isfinite(num) else fallback


def rounded(value, digits=6):
    num = finite_number(value)
    return None if num is None else round(num, digits)


def bool_series(series: pd.Series) -> np.ndarray:
    if series is None:
        return np.zeros(0, dtype=bool)
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).to_numpy(dtype=float) != 0
    return series.fillna("").astype(str).str.lower().isin(["true", "1", "yes", "y"]).to_numpy()


def numeric_series(df: pd.DataFrame, column: str, default: float) -> pd.Series:
    if column not in df:
        return pd.Series(default, index=df.index)
    return pd.to_numeric(df[column], errors="coerce").fillna(default)


def discover_drugs(files: list[Path], compound_types: dict[str, str]) -> list[str]:
    seen: set[str] = set()
    drugs: list[str] = []
    for path in files:
        header = pd.read_csv(path, nrows=0).columns
        for column in header:
            if not column.endswith("_r"):
                continue
            drug = column[:-2]
            if drug.startswith("SRPKIN") or drug not in compound_types:
                continue
            if drug not in seen:
                seen.add(drug)
                drugs.append(drug)
    return drugs


def write_compound_parts(dataset_key: str, out_dir: Path, drug: str, rows: list[list]) -> list[str]:
    paths: list[str] = []
    for part, start in enumerate(range(0, len(rows), PART_SIZE)):
        chunk = rows[start : start + PART_SIZE]
        rel = f"{dataset_key}/compounds/{drug}_part_{part}.json"
        write_json(out_dir / rel, {"drug": drug, "part": part, "rows": chunk})
        paths.append(rel)
    if not paths:
        rel = f"{dataset_key}/compounds/{drug}_part_0.json"
        write_json(out_dir / rel, {"drug": drug, "part": 0, "rows": []})
        paths.append(rel)
    return paths


def write_hover(dataset_key: str, out_dir: Path, drug: str, rows: dict[str, list[float]]) -> str:
    rel = f"{dataset_key}/hover/{drug}.json"
    write_json(out_dir / rel, {"drug": drug, "rows": rows})
    return rel


def make_gene_filename(gene: str) -> str:
    return f"{quote(gene, safe='')}.json"


def clear_dataset_dir(out_dir: Path, dataset_key: str) -> None:
    target = out_dir / dataset_key
    if target.exists():
        shutil.rmtree(target)
    (target / "compounds").mkdir(parents=True, exist_ok=True)
    (target / "hover").mkdir(parents=True, exist_ok=True)
    (target / "sites").mkdir(parents=True, exist_ok=True)


def build_dataset(
    dataset_key: str,
    source_dir: Path,
    out_dir: Path,
    manifest: dict,
    existing_meta: dict[tuple[str, ...], dict],
) -> dict:
    files = sorted(source_dir.glob("*_processed.csv"))
    if not files:
        raise FileNotFoundError(f"No processed CSV files found in {source_dir}")

    compound_types = manifest.get("compoundTypes", {})
    drugs = discover_drugs(files, compound_types)
    drug_set = set(drugs)
    clear_dataset_dir(out_dir, dataset_key)

    row_ids: dict[tuple[str, str, str, str], int] = {}
    catalog: list[dict] = []
    hits_by_row: dict[int, list[list]] = defaultdict(list)
    max_r_by_row: dict[int, float] = defaultdict(lambda: -math.inf)
    hit_count_by_row: dict[int, int] = defaultdict(int)
    observed_count_by_row: dict[int, int] = defaultdict(int)
    drug_hit_counts: dict[str, int] = {drug: 0 for drug in drugs}
    compound_parts: dict[str, list[str]] = {}
    raw_hover_parts: dict[str, str] = {}

    for index, path in enumerate(files, 1):
        header = pd.read_csv(path, nrows=0).columns
        file_drugs = [column[:-2] for column in header if column.endswith("_r") and column[:-2] in drug_set]
        usecols = list(META_COLUMNS)
        default_cols = [column for column in header if column.startswith("default~") and column.endswith("_sn_sum")]
        usecols.extend(default_cols)
        for drug in file_drugs:
            usecols.extend(
                [
                    f"{drug}_r",
                    f"{drug}_pvalue",
                    f"{drug}_mean_sn_cpd",
                    f"{drug}_F_pvalue_dmso",
                    f"{drug}_F_pvalue_cpd",
                    f"{drug}_contamination",
                    f"{drug}_dmso_outlier",
                ]
            )
        df = pd.read_csv(path, usecols=[column for column in usecols if column in header], low_memory=False)
        print(f"[{dataset_key}] {index}/{len(files)} {path.name}: {len(df):,} rows, {len(file_drugs)} compounds")

        current_row_ids: list[int] = []
        for _, source_row in df[META_COLUMNS].iterrows():
            key = source_row_key(source_row)
            row_id = row_ids.get(key)
            if row_id is None:
                row_id = len(catalog)
                row_ids[key] = row_id
                meta = source_metadata(source_row, existing_meta)
                catalog.append(
                    {
                        "i": row_id,
                        **{k: v for k, v in meta.items() if k not in {"x", "y"}},
                        "maxR": 0,
                        "promiscuity": 0,
                        "x": meta.get("x"),
                        "y": meta.get("y"),
                    }
                )
            current_row_ids.append(row_id)

        non_d_cols = [column for column in default_cols if "D_sn_sum" not in column]
        d_cols = [column for column in default_cols if "D_sn_sum" in column]
        dmso_1_series = df[non_d_cols].median(axis=1, skipna=True).fillna(0) if non_d_cols else pd.Series(0, index=df.index)
        dmso_2_series = df[d_cols].median(axis=1, skipna=True).fillna(dmso_1_series) if d_cols else dmso_1_series
        dmso_1 = dmso_1_series.to_numpy(dtype=float)
        dmso_2 = dmso_2_series.to_numpy(dtype=float)

        for drug in file_drugs:
            log2_r = pd.to_numeric(df[f"{drug}_r"], errors="coerce").to_numpy(dtype=float)
            ratio = np.power(2.0, log2_r, where=np.isfinite(log2_r), out=np.full(len(df), np.nan))
            p = pd.to_numeric(df[f"{drug}_pvalue"], errors="coerce").clip(lower=PVALUE_FLOOR, upper=1).to_numpy(dtype=float)
            mean_sn = numeric_series(df, f"{drug}_mean_sn_cpd", 0).to_numpy(dtype=float)
            f_dmso = numeric_series(df, f"{drug}_F_pvalue_dmso", 1).to_numpy(dtype=float)
            f_cpd = numeric_series(df, f"{drug}_F_pvalue_cpd", 1).to_numpy(dtype=float)
            contamination = bool_series(df.get(f"{drug}_contamination"))
            dmso_outlier = bool_series(df.get(f"{drug}_dmso_outlier"))
            if contamination.size == 0:
                contamination = np.zeros(len(df), dtype=bool)
            if dmso_outlier.size == 0:
                dmso_outlier = np.zeros(len(df), dtype=bool)
            hide_dmso = dmso_outlier | (f_dmso < 0.01)
            hide_cpd = contamination | (f_cpd < 0.01)
            finite = np.isfinite(log2_r) & np.isfinite(ratio)

            rows: list[list] = []
            hover: dict[str, list[float]] = {}
            hit_mask = finite & (ratio > HIT_R)
            drug_hit_counts[drug] = int(hit_mask.sum())

            for idx in np.where(finite)[0]:
                row_id = current_row_ids[int(idx)]
                observed_count_by_row[row_id] += 1
                pval = float(p[idx]) if math.isfinite(float(p[idx])) else 1.0
                neglog = -math.log10(max(pval, PVALUE_FLOOR))
                compound_row = [
                    row_id,
                    rounded(ratio[idx]),
                    rounded(pval, 12),
                    rounded(log2_r[idx], 6),
                    rounded(neglog, 6),
                    rounded(mean_sn[idx]),
                    bool(hide_dmso[idx]),
                    bool(hide_cpd[idx]),
                ]
                rows.append(compound_row)
                hover[str(row_id)] = [
                    rounded(dmso_1[idx]),
                    rounded(dmso_2[idx]),
                    rounded(mean_sn[idx]),
                    rounded(mean_sn[idx]),
                ]
                if ratio[idx] > HIT_R:
                    hit = [
                        drug,
                        rounded(ratio[idx]),
                        rounded(pval, 12),
                        rounded(neglog, 6),
                        rounded(pval, 12),
                        rounded(mean_sn[idx]),
                        bool(hide_dmso[idx]),
                        bool(hide_cpd[idx]),
                        compound_types.get(drug, ""),
                    ]
                    hits_by_row[row_id].append(hit)
                    max_r_by_row[row_id] = max(max_r_by_row[row_id], float(ratio[idx]))
                    hit_count_by_row[row_id] += 1

            rows.sort(key=lambda item: item[0])
            compound_parts[drug] = write_compound_parts(dataset_key, out_dir, drug, rows)
            raw_hover_parts[drug] = write_hover(dataset_key, out_dir, drug, hover)

    total_drugs = max(len(drugs), 1)
    for row in catalog:
        row_id = row["i"]
        row_hits = hits_by_row.get(row_id, [])
        row_hits.sort(key=lambda item: item[1], reverse=True)
        max_r = max_r_by_row.get(row_id, 0)
        row["maxR"] = rounded(max_r if math.isfinite(max_r) and max_r > 0 else 0)
        observed = observed_count_by_row.get(row_id, 0)
        row["promiscuity"] = rounded(100 * hit_count_by_row.get(row_id, 0) / max(observed, 1), 12)
        if row.get("x") is None:
            row["x"] = rounded(math.log2(max(row["maxR"], 1e-6))) if row["maxR"] else 0
        if row.get("y") is None:
            row["y"] = row["promiscuity"]

    sites_by_gene: dict[str, list[dict]] = defaultdict(list)
    for row in catalog:
        row_id = row["i"]
        sites_by_gene[row["gene"]].append(
            {
                "row": row_id,
                "label": row["label"],
                "site": row["site"],
                "gene": row["gene"],
                "uniprot": row["uniprot"],
                "positions": row["positions"],
                "maxR": row["maxR"],
                "hits": hits_by_row.get(row_id, []),
            }
        )

    gene_files = {}
    for gene in sorted(sites_by_gene):
        sites = sorted(sites_by_gene[gene], key=lambda item: (-float(item.get("maxR") or 0), item["label"]))
        rel = f"{dataset_key}/sites/{make_gene_filename(gene)}"
        write_json(out_dir / rel, {"gene": gene, "sites": sites})
        gene_files[gene] = rel

    catalog.sort(key=lambda row: row["i"])
    write_json(out_dir / f"{dataset_key}/catalog.json", catalog)

    drug_promiscuity = {
        drug: rounded(100 * drug_hit_counts.get(drug, 0) / max(len(catalog), 1), 12)
        for drug in drugs
    }
    compound_records = [
        {"Promiscuity": drug_promiscuity[drug], "Type": compound_types.get(drug, "")}
        for drug in drugs
    ]
    drug_choices = {
        drug: f"{drug} ({compound_types.get(drug, '')}, {drug_hit_counts.get(drug, 0)} sites hit)"
        for drug in drugs
    }
    default_drug = manifest["datasets"].get(dataset_key, {}).get("defaultDrug")
    if default_drug not in drug_set or drug_hit_counts.get(default_drug, 0) != min(drug_hit_counts.values() or [0]):
        default_drug = min(drugs, key=lambda drug: (drug_hit_counts.get(drug, 0), drugs.index(drug))) if drugs else ""
    default_gene = manifest["datasets"].get(dataset_key, {}).get("defaultGene")
    if default_gene not in gene_files or not any(site.get("maxR", 0) > HIT_R for site in sites_by_gene.get(default_gene, [])):
        default_gene = max(
            sorted(gene_files),
            key=lambda gene: max((site.get("maxR", 0) for site in sites_by_gene.get(gene, [])), default=0),
            default="",
        )
    existing_no_reps = [
        drug for drug in manifest["datasets"].get(dataset_key, {}).get("noReplicateDrugs", []) if drug in drug_set
    ]

    return {
        "label": SOURCE_DATASETS[dataset_key][1],
        "defaultDrug": default_drug,
        "defaultGene": default_gene,
        "rawDrugs": drugs,
        "drugChoices": drug_choices,
        "geneChoices": sorted(gene_files),
        "drugPromiscuity": drug_promiscuity,
        "sitePromiscuity": [row["promiscuity"] for row in catalog],
        "compoundPromiscuityRecords": compound_records,
        "catalog": f"{dataset_key}/catalog.json",
        "compoundParts": compound_parts,
        "geneFiles": gene_files,
        "drugHitCounts": drug_hit_counts,
        "rawHoverParts": raw_hover_parts,
        "activeGeneFilterIndex": None,
        "noReplicateDrugs": existing_no_reps,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        default=r"C:\Users\ethan\Downloads\ypt_statistics\output\best_solution_run_003",
        type=Path,
        help="Directory containing os_0.7 and frac_0.7 processed plex folders.",
    )
    parser.add_argument("--asset-root", default=Path("assets/data"), type=Path)
    args = parser.parse_args()

    manifest_path = args.asset_root / "manifest.json"
    manifest = read_json(manifest_path)
    existing_meta = build_existing_metadata(args.asset_root, manifest)
    new_manifest = dict(manifest)
    new_manifest["datasets"] = {}
    new_manifest["counts"] = {}

    for dataset_key, (source_subdir, _) in SOURCE_DATASETS.items():
        dataset = build_dataset(
            dataset_key,
            args.source_root / source_subdir,
            args.asset_root,
            manifest,
            existing_meta,
        )
        new_manifest["datasets"][dataset_key] = dataset
        new_manifest["counts"][dataset_key] = {
            "rows": len(read_json(args.asset_root / dataset["catalog"])),
            "drugs": len(dataset["rawDrugs"]),
        }

    old_version = str(new_manifest.get("version", "pages-data-v0"))
    match = re.search(r"pages-data-v(\d+)$", old_version)
    if match:
        new_manifest["version"] = f"{old_version[: match.start(1)]}{int(match.group(1)) + 1}"
    else:
        new_manifest["version"] = "pages-data-v12"
    write_json(manifest_path, new_manifest)
    print(f"Wrote {manifest_path} with version {new_manifest['version']}")


if __name__ == "__main__":
    main()
