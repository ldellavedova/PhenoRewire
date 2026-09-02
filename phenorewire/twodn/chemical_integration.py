"""
2D Network annotation expansion for PhenoRewire.

Loads a PhenoRewire rewiring network (GraphML) and a SpecReboot chemical-similarity
network (GraphML), then for each priority feature finds k-hop neighbors in the chemical
network filtered by minimum cosine similarity.

Outputs
-------
TABLE 1 — annotation_summary.csv   : 1 row per priority feature
TABLE 2 — annotation_expansion.csv : 1 row per (feature, neighbor) pair
"""
from __future__ import annotations

import logging
import re
from collections import deque
from pathlib import Path
from typing import Optional

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)

# Column names for table outputs
_SUMMARY_COLS = [
    "feature_id",
    "in_phenorewire",
    "in_specreboot",
    "n_mn_neighbors",
    "n_annotated_neighbors",
    "n_unannotated_neighbors",
    "annotation_confidence",
    "best_neighbor_cosine",
    "top_annotation",
    "rewiring_score",
    "triage_score",
    "hop_depth_used",
    "min_cosine_used",
]

_EXPANSION_COLS = [
    "feature_id",
    "neighbor_id",
    "hop_distance",
    "edge_cosine",
    "neighbor_annotation",
    "neighbor_is_annotated",
    "neighbor_pathway",
    "annotation_confidence",
    "rewiring_score",
    "triage_score",
]

_FIRST_NEIGHBOURS_COLS = [
    "feature_id",
    "rewiring_score",
    "neighbour_id",
    "cosine_similarity",
    "neighbour_annotation",
    "neighbour_pathway",
]


def _load_graph(path) -> nx.Graph:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GraphML not found: {path}")
    G = nx.read_graphml(str(path))
    if G.is_directed():
        G = G.to_undirected()
    return G


def _node_annotation(G: nx.Graph, node: str) -> Optional[str]:
    """Return the annotation label for a node, or None if unannotated."""
    attrs = G.nodes.get(node, {})
    for key in ("consensus_annotation", "annotation", "name", "label", "compound"):
        val = attrs.get(key)
        if val and str(val).strip() and str(val).lower() not in {"nan", "none", "unknown", ""}:
            return str(val).strip()
    return None


def _node_pathway(G: nx.Graph, node: str) -> Optional[str]:
    """Return the NPC pathway/chemical class for a node, or None if absent."""
    attrs = G.nodes.get(node, {})
    for key in ("NPC#pathway", "npc_pathway", "pathway", "chemical_class", "superclass"):
        val = attrs.get(key)
        if val and str(val).strip() and str(val).lower() not in {"nan", "none", "unknown", ""}:
            return str(val).strip()
    return None


def _edge_cosine(G: nx.Graph, u: str, v: str) -> float:
    """Return edge cosine similarity, defaulting to 1.0 if not stored."""
    attrs = G.edges.get((u, v), {})
    for key in ("cosine", "weight", "similarity", "score"):
        val = attrs.get(key)
        if val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                pass
    return 1.0


def _hop_neighbors(
    G: nx.Graph,
    source: str,
    hop_depth: int,
    min_cosine: float,
) -> list[tuple[str, int, float]]:
    """
    BFS up to hop_depth, returning (neighbor_id, hop, min_edge_cosine_on_path).
    Only follows edges with cosine >= min_cosine.
    Excludes source itself.
    """
    if source not in G:
        return []

    visited: dict[str, tuple[int, float]] = {}
    queue: deque[tuple[str, int, float]] = deque([(source, 0, 1.0)])

    while queue:
        node, depth, path_cos = queue.popleft()
        if node in visited:
            continue
        if node != source:
            visited[node] = (depth, path_cos)
        if depth >= hop_depth:
            continue
        for nbr in G.neighbors(node):
            if nbr in visited:
                continue
            cos = _edge_cosine(G, node, nbr)
            if cos < min_cosine:
                continue
            new_cos = min(path_cos, cos)
            queue.append((nbr, depth + 1, new_cos))

    return [(nid, hd, pc) for nid, (hd, pc) in visited.items()]


def build_annotation_expansion(
    phenorewire_graphml: Path,
    specreboot_graphml: Path,
    priority_features_df: pd.DataFrame,
    rewiring_df: pd.DataFrame,
    *,
    top_n: int = 20,
    hop_depth: int = 2,
    min_cosine: float = 0.70,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build annotation expansion tables for priority features.

    Parameters
    ----------
    phenorewire_graphml : path to rewiring_network.graphml produced by PhenoRewire
    specreboot_graphml  : path to SpecReboot molecular network GraphML
    priority_features_df: DataFrame with feature_id (and optionally triage_score/rewiring_score)
                          from priority_features.csv or priority_rewired_nodes.csv
    rewiring_df         : DataFrame with feature_id and rewiring_score; may be the same or empty
    top_n               : maximum priority features to expand
    hop_depth           : maximum BFS depth in SpecReboot network
    min_cosine          : minimum cosine similarity to traverse an edge
    output_dir          : where to write annotation_summary.csv, annotation_expansion.csv,
                          and annotation_first_neighbours.csv

    Returns
    -------
    summary_df, expansion_df, first_neighbours_df
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    G_pr = _load_graph(phenorewire_graphml)
    G_mn = _load_graph(specreboot_graphml)

    logger.info(
        "PhenoRewire network: %d nodes, %d edges",
        G_pr.number_of_nodes(), G_pr.number_of_edges(),
    )
    logger.info(
        "SpecReboot network: %d nodes, %d edges",
        G_mn.number_of_nodes(), G_mn.number_of_edges(),
    )

    # Collect priority feature IDs
    feat_col = _detect_id_col(priority_features_df)
    if feat_col and not priority_features_df.empty:
        feature_ids = priority_features_df[feat_col].astype(str).tolist()[:top_n]
    else:
        feature_ids = list(G_pr.nodes)[:top_n]
        logger.warning("priority_features_df is empty; using first %d PhenoRewire nodes.", top_n)

    # Build lookup dicts for scores
    rew_lookup: dict[str, float] = {}
    rew_id_col = _detect_id_col(rewiring_df)
    if rew_id_col and "rewiring_score" in rewiring_df.columns:
        rew_lookup = dict(zip(
            rewiring_df[rew_id_col].astype(str),
            rewiring_df["rewiring_score"].astype(float),
        ))

    triage_lookup: dict[str, float] = {}
    if feat_col and "triage_score" in priority_features_df.columns:
        triage_lookup = dict(zip(
            priority_features_df[feat_col].astype(str),
            priority_features_df["triage_score"].astype(float),
        ))

    summary_rows = []
    expansion_rows = []
    neighbourhoods: dict[str, list[tuple[str, int, float]]] = {}

    for fid in feature_ids:
        in_pr = fid in G_pr
        in_mn = fid in G_mn
        neighbors = _hop_neighbors(G_mn, fid, hop_depth=hop_depth, min_cosine=min_cosine)
        neighbourhoods[fid] = neighbors

        ann_nbrs = [(nid, h, c) for nid, h, c in neighbors if _node_annotation(G_mn, nid)]
        unann_nbrs = [(nid, h, c) for nid, h, c in neighbors if not _node_annotation(G_mn, nid)]
        n_total = len(neighbors)
        n_ann = len(ann_nbrs)

        annotation_confidence = round(n_ann / n_total, 4) if n_total > 0 else 0.0
        best_cosine = max((c for _, _, c in neighbors), default=float("nan"))
        top_ann = None
        if ann_nbrs:
            best_ann_nbr = max(ann_nbrs, key=lambda x: (x[2], -x[1]))
            top_ann = _node_annotation(G_mn, best_ann_nbr[0])

        summary_rows.append({
            "feature_id": fid,
            "in_phenorewire": in_pr,
            "in_specreboot": in_mn,
            "n_mn_neighbors": n_total,
            "n_annotated_neighbors": n_ann,
            "n_unannotated_neighbors": len(unann_nbrs),
            "annotation_confidence": annotation_confidence,
            "best_neighbor_cosine": round(best_cosine, 4) if not pd.isna(best_cosine) else None,
            "top_annotation": top_ann,
            "rewiring_score": rew_lookup.get(fid, None),
            "triage_score": triage_lookup.get(fid, None),
            "hop_depth_used": hop_depth,
            "min_cosine_used": min_cosine,
        })

        for nid, hop, cos in neighbors:
            ann = _node_annotation(G_mn, nid)
            pathway = _node_pathway(G_mn, nid)
            expansion_rows.append({
                "feature_id": fid,
                "neighbor_id": nid,
                "hop_distance": hop,
                "edge_cosine": round(cos, 4),
                "neighbor_annotation": ann,
                "neighbor_is_annotated": ann is not None,
                "neighbor_pathway": pathway,
                "annotation_confidence": annotation_confidence,
                "rewiring_score": rew_lookup.get(fid, None),
                "triage_score": triage_lookup.get(fid, None),
            })

    summary_df = pd.DataFrame(summary_rows, columns=_SUMMARY_COLS)
    expansion_df = pd.DataFrame(expansion_rows, columns=_EXPANSION_COLS)

    # Write outputs
    summary_df.to_csv(output_dir / "annotation_summary.csv", index=False)
    expansion_df.to_csv(output_dir / "annotation_expansion.csv", index=False)

    # First-neighbours table: hop_distance == 1 only, renamed columns for quick lookup
    first_nbr_df = (
        expansion_df[expansion_df["hop_distance"] == 1]
        .rename(columns={
            "neighbor_id": "neighbour_id",
            "edge_cosine": "cosine_similarity",
            "neighbor_annotation": "neighbour_annotation",
            "neighbor_pathway": "neighbour_pathway",
        })[_FIRST_NEIGHBOURS_COLS]
        .sort_values(["feature_id", "cosine_similarity"], ascending=[True, False])
        .reset_index(drop=True)
    )
    first_nbr_df.to_csv(output_dir / "annotation_first_neighbours.csv", index=False)

    logger.info(
        "Wrote annotation_summary.csv (%d rows), annotation_expansion.csv (%d rows), "
        "annotation_first_neighbours.csv (%d rows) to %s",
        len(summary_df), len(expansion_df), len(first_nbr_df), output_dir,
    )

    # Per-feature subnetwork GraphMLs, reusing the neighbourhoods computed above.
    _write_subnetworks(
        G_mn=G_mn,
        neighbourhoods=neighbourhoods,
        output_dir=output_dir / "subnetworks",
    )

    return summary_df, expansion_df, first_nbr_df


def _detect_id_col(df: pd.DataFrame) -> Optional[str]:
    """Find the feature ID column in a DataFrame."""
    if df is None or df.empty:
        return None
    for candidate in ("feature_id", "node_id", "id", "feature"):
        if candidate in df.columns:
            return candidate
    if len(df.columns) > 0:
        return df.columns[0]
    return None


_UNSAFE_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]|\s')


def _safe_filename(name: str, *, used: set[str]) -> str:
    """Turn a feature ID into a filename that is valid on Windows and unique.

    Feature IDs come from user data and routinely contain characters Windows
    rejects in paths.  Sanitizing alone can also make two distinct IDs collide,
    so a numeric suffix is appended when that happens.
    """
    stem = _UNSAFE_FILENAME_CHARS.sub("_", str(name)).strip(". ") or "feature"
    stem = stem[:100]
    candidate = stem
    suffix = 1
    while candidate.lower() in used:
        suffix += 1
        candidate = f"{stem}_{suffix}"
    used.add(candidate.lower())
    return candidate


def _write_subnetworks(
    G_mn: nx.Graph,
    neighbourhoods: dict[str, list[tuple[str, int, float]]],
    output_dir: Path,
) -> None:
    """Write one GraphML per priority feature containing its k-hop subnetwork."""
    output_dir.mkdir(parents=True, exist_ok=True)
    used: set[str] = set()
    for fid, nbrs in neighbourhoods.items():
        node_set = {fid} | {nid for nid, _, _ in nbrs}
        subG = G_mn.subgraph(node_set).copy()
        out_path = output_dir / f"{_safe_filename(fid, used=used)}_subnetwork.graphml"
        nx.write_graphml(subG, str(out_path))
