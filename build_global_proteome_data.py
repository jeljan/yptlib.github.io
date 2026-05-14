import argparse
import colorsys
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import colormaps
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


APP_DIR = Path(__file__).parent
ASSET_DIR = APP_DIR / "assets" / "data"
DEFAULT_CHECKPOINT_DIR = Path.home() / "Downloads" / "esm3_checkpoints"
DEFAULT_THRESHOLD = 2
DEFAULT_CLUSTER_COLOR = "#8ea9c1"


def as_jsonable(value):
    if value is None:
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        if not math.isfinite(float(value)):
            return None
        return round(float(value), 6)
    return value


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=True, separators=(",", ":"), default=as_jsonable)


def parse_uniprot(value: object) -> str:
    text = str(value or "")
    if "|" in text:
        parts = text.split("|")
        if len(parts) > 1:
            text = parts[1]
    return text.split("-")[0].strip()


def clean_description(value: object) -> str:
    return str(value or "").split(" OS=")[0].strip()


def read_embedding_table(path: Path, source: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing embedding file: {path}")
    df = pd.read_parquet(path)
    df = df.copy()
    df.index = df.index.astype(str)
    df["uniprot"] = df.index.map(parse_uniprot)
    df.index.name = None
    df["embeddingSource"] = source
    df = df[df["uniprot"] != ""]
    return df


def load_embedding_points(checkpoint_dir: Path) -> pd.DataFrame:
    assay = read_embedding_table(checkpoint_dir / "protein_embeddings.parquet", "assay")
    background = read_embedding_table(checkpoint_dir / "human_proteome_embeddings.parquet", "background")
    feature_cols = [col for col in assay.columns if col.startswith("esm3_prot_")]
    combined = pd.concat([assay, background], axis=0, ignore_index=False)
    source_sets = combined.groupby("uniprot")["embeddingSource"].apply(lambda values: sorted(set(values))).to_dict()
    combined["_priority"] = combined["embeddingSource"].map({"assay": 0, "background": 1}).fillna(2)
    combined = combined.sort_values(["uniprot", "_priority"]).drop_duplicates("uniprot", keep="first")
    combined = combined.sort_values(["_priority", "uniprot"]).reset_index(drop=True)
    combined["sourceSet"] = combined["uniprot"].map(source_sets)
    return combined[["uniprot", "embeddingSource", "sourceSet", *feature_cols]]


def load_catalog_overlay(asset_dir: Path) -> dict[str, dict]:
    overlay: dict[str, dict] = defaultdict(
        lambda: {
            "gene": "",
            "description": "",
            "maxRByDataset": {"os": None, "frac": None},
            "hitCountByDataset": {"os": 0, "frac": 0},
            "maxPromiscuityByDataset": {"os": None, "frac": None},
            "contactTypes": set(),
        }
    )
    for dataset_key in ("os", "frac"):
        catalog_path = asset_dir / dataset_key / "catalog.json"
        if not catalog_path.exists():
            continue
        rows = json.loads(catalog_path.read_text(encoding="utf-8"))
        for row in rows:
            uniprot = parse_uniprot(row.get("uniprot") or row.get("proteinId"))
            if not uniprot:
                continue
            entry = overlay[uniprot]
            if not entry["gene"] and row.get("gene"):
                entry["gene"] = str(row["gene"])
            if not entry["description"] and row.get("description"):
                entry["description"] = clean_description(row["description"])
            max_r = row.get("maxR")
            if isinstance(max_r, (int, float)) and math.isfinite(float(max_r)):
                current = entry["maxRByDataset"][dataset_key]
                entry["maxRByDataset"][dataset_key] = max(float(max_r), current or 0.0)
                if float(max_r) >= DEFAULT_THRESHOLD:
                    entry["hitCountByDataset"][dataset_key] += 1
            promiscuity = row.get("promiscuity")
            if isinstance(promiscuity, (int, float)) and math.isfinite(float(promiscuity)):
                current = entry["maxPromiscuityByDataset"][dataset_key]
                entry["maxPromiscuityByDataset"][dataset_key] = max(float(promiscuity), current or 0.0)
            for contact_type in row.get("contactTypes") or []:
                entry["contactTypes"].add(str(contact_type))
    packed = {}
    for uniprot, entry in overlay.items():
        entry["contactTypes"] = sorted(entry["contactTypes"])
        packed[uniprot] = entry
    return packed


def is_seen_assay(assay: dict) -> bool:
    values = (assay.get("maxRByDataset") or {}).values()
    return any(isinstance(value, (int, float)) and math.isfinite(float(value)) for value in values)


def load_domain_metadata(checkpoint_dir: Path) -> dict[str, dict]:
    metadata: dict[str, dict] = {}
    for name in ("domain_log.csv", "human_proteome_domain_log.csv"):
        path = checkpoint_dir / name
        if not path.exists():
            continue
        df = pd.read_csv(path)
        for _, row in df.iterrows():
            uniprot = parse_uniprot(row.get("uniprot"))
            if not uniprot:
                continue
            metadata[uniprot] = {
                "method": str(row.get("method", "") or ""),
                "chunks": int(row.get("n_chunks", 0) or 0),
            }
    for uniprot, annotation in load_interpro_annotations(checkpoint_dir).items():
        metadata.setdefault(uniprot, {})["annotation"] = annotation.get("primary")
        metadata.setdefault(uniprot, {})["annotations"] = annotation.get("all") or []
    return metadata


def optional_download_file(checkpoint_dir: Path, name: str) -> Path | None:
    for candidate in (checkpoint_dir / name, checkpoint_dir.parent / name, Path.home() / "Downloads" / name):
        if candidate.exists():
            return candidate
    return None


def collapse_interpro_labels(df: pd.DataFrame, group_col: str, term_col: str) -> pd.DataFrame:
    clean = df[[group_col, term_col]].dropna().copy()
    roots = clean[term_col].astype(str).str.split(",").str[0]
    roots = roots.str.replace(r"(?i)\s*\b(superfamily|family|domain(s)?|motif|repeat|conserved site)\b\s*", "", regex=True).str.strip()
    clean["root_term"] = roots

    def absorb_substrings(terms: pd.Series) -> str:
        unique_terms = sorted({term for term in terms.dropna().astype(str) if term and term != "nan"}, key=len)
        retained = []
        for term in unique_terms:
            if not any(existing.lower() in term.lower() for existing in retained):
                retained.append(term)
        return "; ".join(retained)

    return clean.groupby(group_col)["root_term"].apply(absorb_substrings).reset_index(name="collapsed_term")


def load_interpro_annotations(checkpoint_dir: Path) -> dict[str, dict]:
    index_path = optional_download_file(checkpoint_dir, "interpro_index.tsv")
    annotation_path = optional_download_file(checkpoint_dir, "annotation_df_interpro_human.tsv")
    if not index_path or not annotation_path:
        return {}

    interpro_index = pd.read_csv(index_path, sep="\t")
    annotation_df = pd.read_csv(annotation_path, sep="\t")
    term_names = dict(zip(interpro_index["ENTRY_AC"].astype(str), interpro_index["ENTRY_NAME"].astype(str)))
    annotation_df = annotation_df.copy()
    annotation_df["uniprot"] = annotation_df["protein_id"].map(parse_uniprot)
    annotation_df["term_name"] = annotation_df["term"].astype(str).map(term_names).fillna(annotation_df["term"].astype(str))
    collapsed = collapse_interpro_labels(annotation_df, "uniprot", "term_name")
    exploded = collapsed.copy()
    exploded["single_term"] = exploded["collapsed_term"].fillna("").str.split("; ")
    exploded = exploded.explode("single_term")
    exploded["single_term"] = exploded["single_term"].astype(str).str.strip()
    exploded = exploded[exploded["single_term"] != ""]
    if exploded.empty:
        return {}

    term_counts = exploded["single_term"].value_counts()
    n_proteins = max(1, collapsed["uniprot"].nunique())
    exploded["idf_score"] = exploded["single_term"].map(lambda term: math.log(n_proteins / term_counts[term]))
    primary = exploded.sort_values(["idf_score", "single_term"], ascending=[False, True]).drop_duplicates("uniprot")
    primary_terms = dict(zip(primary["uniprot"], primary["single_term"]))
    all_terms = {
        row["uniprot"]: [term.strip() for term in str(row["collapsed_term"] or "").split("; ") if term.strip()]
        for _, row in collapsed.iterrows()
    }
    return {
        uniprot: {
            "primary": primary_terms.get(uniprot),
            "all": terms,
        }
        for uniprot, terms in all_terms.items()
    }


def reduce_and_cluster(matrix: np.ndarray, random_state: int):
    scaled = StandardScaler().fit_transform(matrix)
    try:
        from umap import UMAP
        from hdbscan import HDBSCAN

        high_components = min(100, scaled.shape[1], max(2, scaled.shape[0] - 1))
        manifold = UMAP(
            n_components=high_components,
            n_neighbors=50,
            min_dist=0.0,
            metric="cosine",
            random_state=random_state,
        ).fit_transform(scaled)
        labels = HDBSCAN(cluster_selection_epsilon=0.1).fit_predict(manifold)
        coords = UMAP(
            n_components=2,
            n_neighbors=30,
            min_dist=1.0,
            metric="cosine",
            random_state=random_state,
        ).fit_transform(manifold)
        method = "umap-hdbscan"
    except Exception:
        pca_components = min(50, scaled.shape[1], max(2, scaled.shape[0] - 1))
        pca = PCA(n_components=pca_components, svd_solver="randomized", random_state=random_state)
        manifold = pca.fit_transform(scaled)
        cluster_count = max(8, min(80, int(round(math.sqrt(scaled.shape[0] / 4)))))
        labels = MiniBatchKMeans(
            n_clusters=cluster_count,
            random_state=random_state,
            batch_size=8192,
            n_init=5,
        ).fit_predict(manifold)
        coords = manifold[:, :2]
        method = "pca-minibatch-kmeans"
    return coords, labels.astype(str), method


def cycling_spectral_colors(n_clusters: int, n_hue_steps: int = 10) -> list[str]:
    if n_clusters <= 0:
        return []
    n_cycles = int(math.ceil(n_clusters / n_hue_steps))
    cmap = colormaps["Spectral"]
    base_positions = np.linspace(0, 0.8, n_hue_steps, endpoint=True)
    base_rgb = [cmap(float(pos))[:3] for pos in base_positions]
    colors = []
    for cycle in range(n_cycles):
        for r, g, b in base_rgb:
            h, lightness, saturation = colorsys.rgb_to_hls(float(r), float(g), float(b))
            if cycle % 2 == 0:
                new_l = lightness * (1 - 0.3 * (cycle // 2 + 1) / max(n_cycles, 1))
            else:
                new_l = lightness + (1 - lightness) * 0.3 * ((cycle // 2 + 1) / max(n_cycles, 1))
            new_s = np.clip(saturation * (1 - 0.15 * cycle / max(n_cycles - 1, 1)), 0.3, 1.0)
            nr, ng, nb = colorsys.hls_to_rgb(h, float(np.clip(new_l, 0.15, 0.9)), float(new_s))
            colors.append(f"rgb({int(nr * 255)},{int(ng * 255)},{int(nb * 255)})")
            if len(colors) >= n_clusters:
                return colors
    return colors[:n_clusters]


def sort_annotations_by_cluster_prevalence(points: list[dict], cluster_counts: dict[str, int]) -> dict[str, dict]:
    term_counts_by_cluster: dict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
    for point in points:
        cluster = str(point.get("cluster", "-1"))
        terms = {str(term).strip() for term in point.get("annotations") or [] if str(term).strip()}
        for term in terms:
            term_counts_by_cluster[cluster][term] += 1

    for point in points:
        cluster = str(point.get("cluster", "-1"))
        cluster_size = max(1, int(cluster_counts.get(cluster, 0)))
        terms = [str(term).strip() for term in point.get("annotations") or [] if str(term).strip()]
        terms = sorted(
            dict.fromkeys(terms),
            key=lambda term: (-(term_counts_by_cluster[cluster][term] / cluster_size), term.lower()),
        )
        point["annotations"] = terms
        point["annotation"] = terms[0] if terms else None

    summaries = {}
    for cluster, term_counts in term_counts_by_cluster.items():
        cluster_size = max(1, int(cluster_counts.get(cluster, 0)))
        if not term_counts:
            continue
        top_terms = sorted(term_counts.items(), key=lambda item: (-(item[1] / cluster_size), item[0].lower()))[:2]
        term, count = top_terms[0]
        summaries[cluster] = {
            "annotation": term,
            "annotationPercent": round((count / cluster_size) * 100, 1),
            "annotationCount": count,
            "topAnnotations": [
                {
                    "annotation": top_term,
                    "percent": round((top_count / cluster_size) * 100, 1),
                    "count": top_count,
                }
                for top_term, top_count in top_terms
            ],
        }
    return summaries


def update_manifest(asset_dir: Path, relative_path: str):
    manifest_path = asset_dir / "manifest.json"
    if not manifest_path.exists():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["globalProteome"] = relative_path
    write_json(manifest_path, manifest)


def build_payload(checkpoint_dir: Path, asset_dir: Path, random_state: int):
    embeddings = load_embedding_points(checkpoint_dir)
    overlay = load_catalog_overlay(asset_dir)
    domains = load_domain_metadata(checkpoint_dir)
    feature_cols = [col for col in embeddings.columns if col.startswith("esm3_prot_")]
    matrix = embeddings[feature_cols].to_numpy(dtype=np.float32, copy=False)
    coords, labels, method = reduce_and_cluster(matrix, random_state)

    cluster_ids = sorted(set(labels), key=lambda value: (value == "-1", int(value) if value.lstrip("-").isdigit() else value))
    non_noise = [cluster_id for cluster_id in cluster_ids if cluster_id != "-1"]
    color_map = {"-1": "rgba(200,200,200,0.3)"}
    color_map.update({cluster_id: color for cluster_id, color in zip(non_noise, cycling_spectral_colors(len(non_noise)))})

    points = []
    cluster_counts = defaultdict(int)
    for idx, row in embeddings.iterrows():
        uniprot = row["uniprot"]
        assay = overlay.get(uniprot, {})
        domain = domains.get(uniprot, {})
        cluster = labels[idx]
        cluster_counts[cluster] += 1
        source_parts = set(row.get("sourceSet") or [])
        if uniprot in overlay:
            source_parts.add("assay")
        source = "both" if source_parts == {"assay", "background"} else next(iter(source_parts or {"background"}))
        annotation = domain.get("annotation")
        if not annotation and domain.get("method"):
            annotation = f"{domain['method']} domains" if domain.get("chunks", 0) else domain["method"]
        annotations = domain.get("annotations") or ([annotation] if annotation else [])
        points.append(
            {
                "uniprot": uniprot,
                "gene": assay.get("gene") or uniprot,
                "x": float(coords[idx, 0]),
                "y": float(coords[idx, 1]),
                "cluster": cluster,
                "annotation": annotation,
                "annotations": annotations,
                "source": source,
                "seen": is_seen_assay(assay),
                "description": assay.get("description") or "",
                "maxRByDataset": assay.get("maxRByDataset") or {"os": None, "frac": None},
                "hitCountByDataset": assay.get("hitCountByDataset") or {"os": 0, "frac": 0},
                "maxPromiscuityByDataset": assay.get("maxPromiscuityByDataset") or {"os": None, "frac": None},
                "contactTypes": assay.get("contactTypes") or [],
            }
        )

    cluster_annotations = sort_annotations_by_cluster_prevalence(points, cluster_counts)

    clusters = [
        {
            "id": cluster_id,
            "label": "Noise" if cluster_id == "-1" else f"Cluster {cluster_id}",
            "color": color_map.get(cluster_id, "rgba(200,200,200,0.3)"),
            "count": cluster_counts[cluster_id],
            **cluster_annotations.get(cluster_id, {}),
        }
        for cluster_id in cluster_ids
    ]
    return {
        "version": 1,
        "method": method,
        "points": points,
        "clusters": clusters,
        "defaults": {"threshold": DEFAULT_THRESHOLD},
    }


def rebuild_payload_with_existing_layout(checkpoint_dir: Path, asset_dir: Path, layout_path: Path):
    existing = json.loads(layout_path.read_text(encoding="utf-8"))
    layout_by_uniprot = {parse_uniprot(point.get("uniprot")): point for point in existing.get("points", [])}
    overlay = load_catalog_overlay(asset_dir)
    domains = load_domain_metadata(checkpoint_dir)
    points = []
    cluster_counts = defaultdict(int)

    for uniprot, layout in layout_by_uniprot.items():
        if not uniprot:
            continue
        assay = overlay.get(uniprot, {})
        domain = domains.get(uniprot, {})
        annotation = domain.get("annotation")
        if not annotation and domain.get("method"):
            annotation = f"{domain['method']} domains" if domain.get("chunks", 0) else domain["method"]
        annotations = domain.get("annotations") or ([annotation] if annotation else [])
        cluster = str(layout.get("cluster", "-1"))
        cluster_counts[cluster] += 1
        points.append(
            {
                "uniprot": uniprot,
                "gene": assay.get("gene") or layout.get("gene") or uniprot,
                "x": float(layout["x"]),
                "y": float(layout["y"]),
                "cluster": cluster,
                "annotation": annotation,
                "annotations": annotations,
                "source": layout.get("source") or ("assay" if uniprot in overlay else "background"),
                "seen": is_seen_assay(assay),
                "description": assay.get("description") or layout.get("description") or "",
                "maxRByDataset": assay.get("maxRByDataset") or {"os": None, "frac": None},
                "hitCountByDataset": assay.get("hitCountByDataset") or {"os": 0, "frac": 0},
                "maxPromiscuityByDataset": assay.get("maxPromiscuityByDataset") or {"os": None, "frac": None},
                "contactTypes": assay.get("contactTypes") or [],
            }
        )

    cluster_annotations = sort_annotations_by_cluster_prevalence(points, cluster_counts)

    cluster_color_by_id = {str(cluster.get("id")): cluster.get("color") for cluster in existing.get("clusters", [])}
    cluster_ids = sorted(cluster_counts, key=lambda value: (value == "-1", int(value) if value.lstrip("-").isdigit() else value))
    clusters = [
        {
            "id": cluster_id,
            "label": "Noise" if cluster_id == "-1" else f"Cluster {cluster_id}",
            "color": cluster_color_by_id.get(cluster_id) or ("rgba(200,200,200,0.3)" if cluster_id == "-1" else DEFAULT_CLUSTER_COLOR),
            "count": cluster_counts[cluster_id],
            **cluster_annotations.get(cluster_id, {}),
        }
        for cluster_id in cluster_ids
    ]
    return {
        "version": 1,
        "method": f"{existing.get('method', 'external-layout')}-metadata-refresh",
        "points": points,
        "clusters": clusters,
        "defaults": {"threshold": DEFAULT_THRESHOLD},
    }


def payload_from_coordinate_csv(checkpoint_dir: Path, asset_dir: Path, coordinate_csv: Path):
    coords_df = pd.read_csv(coordinate_csv)
    required = {"UMAP1", "UMAP2", "cluster", "protein_id"}
    missing = required - set(coords_df.columns)
    if missing:
        raise ValueError(f"Coordinate CSV is missing required columns: {', '.join(sorted(missing))}")
    coords_df = coords_df.copy()
    coords_df["uniprot"] = coords_df["protein_id"].map(parse_uniprot)
    coords_df = coords_df[coords_df["uniprot"] != ""]
    coords_df = coords_df.drop_duplicates("uniprot", keep="first")

    overlay = load_catalog_overlay(asset_dir)
    domains = load_domain_metadata(checkpoint_dir)
    embeddings = load_embedding_points(checkpoint_dir)
    source_sets = {row["uniprot"]: set(row.get("sourceSet") or []) for _, row in embeddings.iterrows()}
    cluster_counts = defaultdict(int)
    points = []

    for _, layout in coords_df.iterrows():
        uniprot = layout["uniprot"]
        assay = overlay.get(uniprot, {})
        domain = domains.get(uniprot, {})
        annotation = domain.get("annotation")
        if not annotation and domain.get("method"):
            annotation = f"{domain['method']} domains" if domain.get("chunks", 0) else domain["method"]
        annotations = domain.get("annotations") or ([annotation] if annotation else [])
        cluster = str(layout.get("cluster", "-1"))
        cluster_counts[cluster] += 1
        source_parts = set(source_sets.get(uniprot) or [])
        if uniprot in overlay:
            source_parts.add("assay")
        source = "both" if source_parts == {"assay", "background"} else next(iter(source_parts or {"background"}))
        points.append(
            {
                "uniprot": uniprot,
                "gene": assay.get("gene") or (str(layout.get("gene") or "").strip() or uniprot),
                "x": float(layout["UMAP1"]),
                "y": float(layout["UMAP2"]),
                "cluster": cluster,
                "annotation": annotation,
                "annotations": annotations,
                "source": source,
                "seen": is_seen_assay(assay),
                "description": assay.get("description") or "",
                "maxRByDataset": assay.get("maxRByDataset") or {"os": None, "frac": None},
                "hitCountByDataset": assay.get("hitCountByDataset") or {"os": 0, "frac": 0},
                "maxPromiscuityByDataset": assay.get("maxPromiscuityByDataset") or {"os": None, "frac": None},
                "contactTypes": assay.get("contactTypes") or [],
            }
        )

    cluster_annotations = sort_annotations_by_cluster_prevalence(points, cluster_counts)

    cluster_ids = sorted(cluster_counts, key=lambda value: (value == "-1", int(value) if value.lstrip("-").isdigit() else value))
    non_noise = [cluster_id for cluster_id in cluster_ids if cluster_id != "-1"]
    color_map = {"-1": "rgba(200,200,200,0.3)"}
    color_map.update({cluster_id: color for cluster_id, color in zip(non_noise, cycling_spectral_colors(len(non_noise), n_hue_steps=7))})
    clusters = [
        {
            "id": cluster_id,
            "label": "Noise" if cluster_id == "-1" else f"Cluster {cluster_id}",
            "color": color_map.get(cluster_id, DEFAULT_CLUSTER_COLOR),
            "count": cluster_counts[cluster_id],
            **cluster_annotations.get(cluster_id, {}),
        }
        for cluster_id in cluster_ids
    ]
    return {
        "version": 1,
        "method": f"external-coordinate-csv:{coordinate_csv.name}",
        "points": points,
        "clusters": clusters,
        "defaults": {"threshold": DEFAULT_THRESHOLD},
    }


def main():
    parser = argparse.ArgumentParser(description="Build the global ESM3 proteome visualization asset.")
    parser.add_argument("--checkpoint-dir", type=Path, default=DEFAULT_CHECKPOINT_DIR)
    parser.add_argument("--asset-dir", type=Path, default=ASSET_DIR)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--coordinate-csv",
        type=Path,
        default=None,
        help="Use precomputed coordinates/clusters from a CSV with UMAP1, UMAP2, cluster, protein_id columns.",
    )
    parser.add_argument(
        "--reuse-existing-layout",
        action="store_true",
        help="Refresh metadata/annotations/overlays while preserving the existing global_proteome.json coordinates and clusters.",
    )
    args = parser.parse_args()

    output_path = args.asset_dir / "global_proteome.json"
    if args.coordinate_csv:
        payload = payload_from_coordinate_csv(args.checkpoint_dir, args.asset_dir, args.coordinate_csv)
    elif args.reuse_existing_layout:
        if not output_path.exists():
            raise FileNotFoundError(f"Cannot reuse layout because {output_path} does not exist")
        payload = rebuild_payload_with_existing_layout(args.checkpoint_dir, args.asset_dir, output_path)
    else:
        payload = build_payload(args.checkpoint_dir, args.asset_dir, args.random_state)
    write_json(output_path, payload)
    update_manifest(args.asset_dir, "global_proteome.json")
    print(
        f"global proteome asset written to {output_path} "
        f"({len(payload['points']):,} points, {len(payload['clusters']):,} clusters, {payload['method']})"
    )


if __name__ == "__main__":
    main()
