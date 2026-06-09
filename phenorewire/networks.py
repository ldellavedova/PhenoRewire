# ── PhenoRewire · Step 3 — Network construction & rewiring ───────────────────
# Builds one Spearman correlation network per biological state on the selected
# feature set. Computes topological metrics: degree, degree centrality,
# weighted degree (abs rho), betweenness centrality (distance = 1/abs(rho)),
# eigenvector centrality, Louvain community membership, network modularity.
# Quantifies edge-by-edge rewiring: shared / state_A_only / state_B_only / sign_switch.
from __future__ import annotations

from pathlib import Path
from typing import Tuple, Optional, Dict, Any

import logging
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

EDGE_COLS = ["node_a", "node_b", "r", "q", "sign"]
NODE_METRIC_COLS = [
    "degree",
    "degree_centrality",
    "weighted_degree_abs",
    "betweenness_centrality",
    "eigenvector_centrality",
    "community_id",
    "community_size",
]


def _threshold_grid(base_abs_r_thr: float, base_fdr_alpha: float, floor_abs_r_thr: float, ceil_fdr_alpha: float) -> list[tuple[float, float]]:
    abs_candidates = sorted(
        {
            round(min(0.95, base_abs_r_thr + 0.10), 4),
            round(min(0.95, base_abs_r_thr + 0.05), 4),
            round(base_abs_r_thr, 4),
            round(max(floor_abs_r_thr, base_abs_r_thr - 0.05), 4),
            round(max(floor_abs_r_thr, base_abs_r_thr - 0.10), 4),
            round(floor_abs_r_thr, 4),
        },
        reverse=True,
    )
    fdr_candidates = sorted(
        {
            round(min(base_fdr_alpha, 0.01), 4),
            round(base_fdr_alpha, 4),
            round(min(ceil_fdr_alpha, max(base_fdr_alpha, 0.10)), 4),
            round(ceil_fdr_alpha, 4),
        }
    )
    return [(float(a), float(q)) for a in abs_candidates for q in fdr_candidates]


def _select_adaptive_thresholds(
    rvals: np.ndarray,
    qvals: np.ndarray,
    *,
    n_nodes: int,
    base_abs_r_thr: float,
    base_fdr_alpha: float,
    adaptive: bool,
    min_edges: int,
    max_density: float,
    floor_abs_r_thr: float,
    ceil_fdr_alpha: float,
) -> dict[str, Any]:
    total_possible_edges = int(n_nodes * (n_nodes - 1) / 2)
    target_min_edges = min(int(min_edges), total_possible_edges) if total_possible_edges > 0 else 0
    target_max_edges = int(np.floor(float(max_density) * total_possible_edges)) if total_possible_edges > 0 else 0
    target_max_edges = max(target_min_edges, target_max_edges)

    finite = np.isfinite(rvals) & np.isfinite(qvals)
    candidate_edges_before_filter = int(finite.sum())
    if candidate_edges_before_filter == 0 or total_possible_edges == 0:
        return {
            "base_abs_r_thr": float(base_abs_r_thr),
            "base_fdr_alpha": float(base_fdr_alpha),
            "applied_abs_r_thr": float(base_abs_r_thr),
            "applied_fdr_alpha": float(base_fdr_alpha),
            "adaptive_thresholds_used": False,
            "target_min_edges": int(target_min_edges),
            "target_max_edges": int(target_max_edges),
            "candidate_edges_before_filter": int(candidate_edges_before_filter),
            "selected_edge_count": 0,
            "selected_density": 0.0,
        }

    if not adaptive:
        mask = finite & (qvals <= float(base_fdr_alpha)) & (np.abs(rvals) >= float(base_abs_r_thr))
        return {
            "base_abs_r_thr": float(base_abs_r_thr),
            "base_fdr_alpha": float(base_fdr_alpha),
            "applied_abs_r_thr": float(base_abs_r_thr),
            "applied_fdr_alpha": float(base_fdr_alpha),
            "adaptive_thresholds_used": False,
            "target_min_edges": int(target_min_edges),
            "target_max_edges": int(target_max_edges),
            "candidate_edges_before_filter": int(candidate_edges_before_filter),
            "selected_edge_count": int(mask.sum()),
            "selected_density": float(mask.sum() / total_possible_edges) if total_possible_edges > 0 else 0.0,
        }

    best_choice = None
    best_score = None
    for abs_thr, fdr_thr in _threshold_grid(base_abs_r_thr, base_fdr_alpha, floor_abs_r_thr, ceil_fdr_alpha):
        mask = finite & (qvals <= float(fdr_thr)) & (np.abs(rvals) >= float(abs_thr))
        edge_count = int(mask.sum())
        density = float(edge_count / total_possible_edges) if total_possible_edges > 0 else 0.0
        too_few = max(0, target_min_edges - edge_count)
        too_many = max(0, edge_count - target_max_edges)
        deviation = (abs(abs_thr - base_abs_r_thr) / 0.05) + (abs(fdr_thr - base_fdr_alpha) / max(base_fdr_alpha, 1e-6))
        score = (too_few * 1000.0) + (too_many * 10.0) + deviation + density
        choice = {
            "base_abs_r_thr": float(base_abs_r_thr),
            "base_fdr_alpha": float(base_fdr_alpha),
            "applied_abs_r_thr": float(abs_thr),
            "applied_fdr_alpha": float(fdr_thr),
            "adaptive_thresholds_used": bool(abs(abs_thr - base_abs_r_thr) > 1e-12 or abs(fdr_thr - base_fdr_alpha) > 1e-12),
            "target_min_edges": int(target_min_edges),
            "target_max_edges": int(target_max_edges),
            "candidate_edges_before_filter": int(candidate_edges_before_filter),
            "selected_edge_count": int(edge_count),
            "selected_density": float(density),
        }
        if best_choice is None or score < best_score:
            best_choice = choice
            best_score = score

    return best_choice


# ============================================================
# Correlation Network
# ============================================================

def correlation_network(
    X_feat_by_samp: pd.DataFrame,
    *,
    abs_r_thr: float,
    fdr_alpha: float,
    adaptive_thresholds: bool = True,
    min_edges: int = 3,
    max_density: float = 0.25,
    floor_abs_r_thr: float = 0.5,
    ceil_fdr_alpha: float = 0.2,
) -> Tuple[nx.Graph, pd.DataFrame, dict[str, Any]]:
    """
    Build an undirected Spearman correlation network among features.

    Parameters
    ----------
    X_feat_by_samp
        DataFrame (features x samples). Index = feature_id; columns = sample_id.
    abs_r_thr
        Absolute Spearman rho threshold applied after FDR filtering.
    fdr_alpha
        BH FDR alpha for edge p-values.

    Returns
    -------
    G, edges_df
        NetworkX graph and an edges table with columns:
        node_a, node_b, r, q, sign
    """

    # Basic shape checks
    if X_feat_by_samp is None or X_feat_by_samp.empty:
        return nx.Graph(), pd.DataFrame(columns=EDGE_COLS), {"adaptive_thresholds_used": False, "applied_abs_r_thr": float(abs_r_thr), "applied_fdr_alpha": float(fdr_alpha), "selected_edge_count": 0, "selected_density": 0.0, "target_min_edges": 0, "target_max_edges": 0, "candidate_edges_before_filter": 0, "base_abs_r_thr": float(abs_r_thr), "base_fdr_alpha": float(fdr_alpha)}

    if X_feat_by_samp.shape[0] < 2 or X_feat_by_samp.shape[1] < 3:
        # need >=2 features and >=3 samples for a meaningful rank correlation matrix
        return nx.Graph(), pd.DataFrame(columns=EDGE_COLS), {"adaptive_thresholds_used": False, "applied_abs_r_thr": float(abs_r_thr), "applied_fdr_alpha": float(fdr_alpha), "selected_edge_count": 0, "selected_density": 0.0, "target_min_edges": 0, "target_max_edges": 0, "candidate_edges_before_filter": 0, "base_abs_r_thr": float(abs_r_thr), "base_fdr_alpha": float(fdr_alpha)}

    # Raw matrix (features x samples)
    mat = X_feat_by_samp.values.astype(float)

    # Replace inf with NaN (do NOT force NaN->0; that can create artificial constants)
    mat[~np.isfinite(mat)] = np.nan

    # Drop constant/invalid features within this group (nanstd ignores NaNs)
    std = np.nanstd(mat, axis=1)
    keep = np.isfinite(std) & (std > 0)

    n_drop = int((~keep).sum())
    if n_drop > 0:
        logger.warning(
            "Dropping %d/%d constant (or invalid) features before correlation.",
            n_drop, mat.shape[0],
        )

    mat2 = mat[keep, :]
    feats2 = X_feat_by_samp.index.astype(str).to_numpy()[keep].tolist()

    # Need enough features & samples after filtering
    if mat2.shape[0] < 2 or mat2.shape[1] < 3:
        return nx.Graph(), pd.DataFrame(columns=EDGE_COLS), {"adaptive_thresholds_used": False, "applied_abs_r_thr": float(abs_r_thr), "applied_fdr_alpha": float(fdr_alpha), "selected_edge_count": 0, "selected_density": 0.0, "target_min_edges": 0, "target_max_edges": 0, "candidate_edges_before_filter": 0, "base_abs_r_thr": float(abs_r_thr), "base_fdr_alpha": float(fdr_alpha)}

    # Memory note: O(n_features²) — ~8 MB at 1k features, ~200 MB at 5k,
    # ~800 MB at 10k. For large feature sets consider pre-filtering.
    R, P = spearmanr(mat2.T, axis=0)

    if not isinstance(R, np.ndarray) or R.ndim != 2:
        return nx.Graph(), pd.DataFrame(columns=EDGE_COLS), {"adaptive_thresholds_used": False, "applied_abs_r_thr": float(abs_r_thr), "applied_fdr_alpha": float(fdr_alpha), "selected_edge_count": 0, "selected_density": 0.0, "target_min_edges": 0, "target_max_edges": 0, "candidate_edges_before_filter": 0, "base_abs_r_thr": float(abs_r_thr), "base_fdr_alpha": float(fdr_alpha)}

    # Upper triangle indices
    iu = np.triu_indices(R.shape[0], k=1)
    rvals = R[iu].astype(float)
    pvals = P[iu].astype(float)

    finite = np.isfinite(rvals) & np.isfinite(pvals)
    if int(finite.sum()) == 0:
        return nx.Graph(), pd.DataFrame(columns=EDGE_COLS), {"adaptive_thresholds_used": False, "applied_abs_r_thr": float(abs_r_thr), "applied_fdr_alpha": float(fdr_alpha), "selected_edge_count": 0, "selected_density": 0.0, "target_min_edges": 0, "target_max_edges": 0, "candidate_edges_before_filter": 0, "base_abs_r_thr": float(abs_r_thr), "base_fdr_alpha": float(fdr_alpha)}

    # BH-FDR on p-values
    rej = np.zeros_like(pvals, dtype=bool)
    qvals = np.full_like(pvals, np.nan, dtype=float)

    rej_f, q_f, *_ = multipletests(
        pvals[finite],
        alpha=float(fdr_alpha),
        method="fdr_bh",
    )
    rej[finite] = rej_f
    qvals[finite] = q_f

    threshold_meta = _select_adaptive_thresholds(
        rvals,
        qvals,
        n_nodes=len(feats2),
        base_abs_r_thr=float(abs_r_thr),
        base_fdr_alpha=float(fdr_alpha),
        adaptive=bool(adaptive_thresholds),
        min_edges=int(min_edges),
        max_density=float(max_density),
        floor_abs_r_thr=float(floor_abs_r_thr),
        ceil_fdr_alpha=float(ceil_fdr_alpha),
    )

    mask = (
        np.isfinite(rvals)
        & np.isfinite(qvals)
        & (qvals <= float(threshold_meta["applied_fdr_alpha"]))
        & (np.abs(rvals) >= float(threshold_meta["applied_abs_r_thr"]))
    )

    a_idx = iu[0][mask]
    b_idx = iu[1][mask]

    if a_idx.size == 0:
        threshold_meta["selected_edge_count"] = 0
        threshold_meta["selected_density"] = 0.0
        return nx.Graph(), pd.DataFrame(columns=EDGE_COLS), threshold_meta

    edges_df = pd.DataFrame(
        {
            "node_a": [feats2[i] for i in a_idx],
            "node_b": [feats2[j] for j in b_idx],
            "r": rvals[mask],
            "q": qvals[mask],
        }
    )
    edges_df["sign"] = np.where(edges_df["r"].astype(float) >= 0, "pos", "neg")

    # Build graph
    G = nx.Graph()
    for a, b, r, q, sign in edges_df.itertuples(index=False):
        G.add_edge(str(a), str(b), r=float(r), q=float(q), sign=str(sign))

    return G, edges_df, threshold_meta


def _build_weighted_graph(edges_df: Optional[pd.DataFrame]) -> nx.Graph:
    G = nx.Graph()
    if edges_df is None or edges_df.empty:
        return G

    for a, b, r, q, sign in edges_df[EDGE_COLS].itertuples(index=False):
        weight_abs = abs(float(r))
        G.add_edge(
            str(a),
            str(b),
            r=float(r),
            q=float(q),
            sign=str(sign),
            weight_abs=weight_abs,
            distance=(1.0 / max(weight_abs, 1e-12)),
        )
    return G


def _eigenvector_by_component(G: nx.Graph) -> Dict[str, float]:
    scores = {str(node): 0.0 for node in G.nodes()}
    if G.number_of_nodes() == 0:
        return scores

    for nodes in nx.connected_components(G):
        sub = G.subgraph(nodes).copy()
        if sub.number_of_nodes() == 1 or sub.number_of_edges() == 0:
            only = next(iter(sub.nodes()))
            scores[str(only)] = 0.0
            continue
        try:
            comp_scores = nx.eigenvector_centrality(sub, weight="weight_abs", max_iter=2000)
        except Exception:
            comp_scores = {str(node): 0.0 for node in sub.nodes()}
        for node, value in comp_scores.items():
            scores[str(node)] = float(value)
    return scores


def _community_detection(G: nx.Graph, seed: int = 42) -> tuple[Dict[str, int], pd.DataFrame, Dict[str, Any]]:
    if G.number_of_nodes() == 0:
        empty = pd.DataFrame(columns=["community_id", "size", "n_edges", "members"])
        meta = {"algorithm": "louvain", "modularity": 0.0, "n_communities": 0}
        return {}, empty, meta

    try:
        communities = nx.community.louvain_communities(G, weight="weight_abs", seed=int(seed))
        algorithm = "louvain"
    except Exception:
        communities = list(nx.community.greedy_modularity_communities(G, weight="weight_abs"))
        algorithm = "greedy_modularity"

    membership: Dict[str, int] = {}
    records = []
    for cid, members in enumerate(communities, start=1):
        member_list = sorted(str(node) for node in members)
        for node in member_list:
            membership[node] = cid
        sub = G.subgraph(member_list)
        records.append(
            {
                "community_id": cid,
                "size": len(member_list),
                "n_edges": int(sub.number_of_edges()),
                "members": ";".join(member_list),
            }
        )

    modularity = 0.0
    if communities and G.number_of_edges() > 0:
        modularity = float(nx.community.modularity(G, communities, weight="weight_abs"))

    summary = pd.DataFrame(records).sort_values(["size", "n_edges"], ascending=[False, False])
    meta = {
        "algorithm": algorithm,
        "modularity": modularity,
        "n_communities": len(communities),
    }
    return membership, summary, meta


def compute_network_metrics(
    selected_features: pd.DataFrame,
    edges_df: Optional[pd.DataFrame],
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    nodes_df = selected_features.copy()
    nodes_df.index = nodes_df.index.astype(str)

    G = _build_weighted_graph(edges_df)
    if not nodes_df.empty:
        G.add_nodes_from(nodes_df.index.tolist())

    if G.number_of_nodes() == 0:
        for col in NODE_METRIC_COLS:
            nodes_df[col] = []
        summary = pd.DataFrame(
            [
                {
                    "n_nodes": 0,
                    "n_edges": 0,
                    "density": 0.0,
                    "n_components": 0,
                    "largest_component_size": 0,
                    "positive_edges": 0,
                    "negative_edges": 0,
                    "mean_abs_r": 0.0,
                    "community_algorithm": "louvain",
                    "n_communities": 0,
                    "modularity": 0.0,
                }
            ]
        )
        return nodes_df, pd.DataFrame(columns=["community_id", "size", "n_edges", "members"]), summary

    degree = dict(G.degree())
    degree_cent = nx.degree_centrality(G) if G.number_of_nodes() > 1 else {n: 0.0 for n in G.nodes()}
    weighted_degree = dict(G.degree(weight="weight_abs"))

    if G.number_of_edges() == 0:
        betweenness = {str(node): 0.0 for node in G.nodes()}
    else:
        betweenness = nx.betweenness_centrality(G, weight="distance", normalized=True)

    eigenvector = _eigenvector_by_component(G)
    membership, community_summary, community_meta = _community_detection(G, seed=seed)
    community_sizes = {
        int(row["community_id"]): int(row["size"])
        for _, row in community_summary.iterrows()
    }

    nodes_df["degree"] = nodes_df.index.map(lambda x: int(degree.get(str(x), 0)))
    nodes_df["degree_centrality"] = nodes_df.index.map(lambda x: float(degree_cent.get(str(x), 0.0)))
    nodes_df["weighted_degree_abs"] = nodes_df.index.map(lambda x: float(weighted_degree.get(str(x), 0.0)))
    nodes_df["betweenness_centrality"] = nodes_df.index.map(lambda x: float(betweenness.get(str(x), 0.0)))
    nodes_df["eigenvector_centrality"] = nodes_df.index.map(lambda x: float(eigenvector.get(str(x), 0.0)))
    nodes_df["community_id"] = nodes_df.index.map(lambda x: int(membership.get(str(x), -1)))
    nodes_df["community_size"] = nodes_df["community_id"].map(lambda cid: int(community_sizes.get(int(cid), 1 if cid != -1 else 0)))

    components = list(nx.connected_components(G))
    largest_component = max((len(comp) for comp in components), default=0)
    mean_abs_r = 0.0
    pos_edges = 0
    neg_edges = 0
    if edges_df is not None and not edges_df.empty:
        mean_abs_r = float(edges_df["r"].abs().mean())
        pos_edges = int((edges_df["sign"].astype(str) == "pos").sum())
        neg_edges = int((edges_df["sign"].astype(str) == "neg").sum())

    summary = pd.DataFrame(
        [
            {
                "n_nodes": int(G.number_of_nodes()),
                "n_edges": int(G.number_of_edges()),
                "density": float(nx.density(G)) if G.number_of_nodes() > 1 else 0.0,
                "n_components": int(nx.number_connected_components(G)),
                "largest_component_size": int(largest_component),
                "positive_edges": pos_edges,
                "negative_edges": neg_edges,
                "mean_abs_r": mean_abs_r,
                "community_algorithm": community_meta["algorithm"],
                "n_communities": int(community_meta["n_communities"]),
                "modularity": float(community_meta["modularity"]),
            }
        ]
    )
    return nodes_df, community_summary, summary


# ============================================================
# Louvain Stability Check
# ============================================================

def check_louvain_stability(
    edges_df: Optional[pd.DataFrame],
    *,
    n_runs: int = 10,
    base_seed: int = 42,
) -> dict[str, Any]:
    """
    Run Louvain community detection ``n_runs`` times with different seeds and
    report the mean adjusted mutual information (AMI) between all pairs of runs.

    Parameters
    ----------
    edges_df
        Edge DataFrame as produced by ``correlation_network()``.
    n_runs
        Number of independent Louvain runs.  10 is sufficient for a stability estimate.
    base_seed
        First seed; subsequent seeds are ``base_seed + i``.

    Returns
    -------
    dict with keys:
        ami_mean, ami_min, n_runs, stable (bool, True if ami_mean >= 0.8),
        warning (str or None)

    Example
    -------
    >>> meta = check_louvain_stability(edges_df, n_runs=10, base_seed=42)
    >>> if not meta['stable']:
    ...     print(meta['warning'])
    """
    try:
        from sklearn.metrics import adjusted_mutual_info_score
    except ImportError as exc:
        raise ImportError(
            "LOUVAIN_STABILITY_CHECK requires scikit-learn. "
            "Install it with: pip install phenorewire[benchmark]"
        ) from exc

    G = _build_weighted_graph(edges_df)
    if G.number_of_nodes() < 2 or G.number_of_edges() == 0:
        return {
            "ami_mean": 1.0,
            "ami_min": 1.0,
            "n_runs": 0,
            "stable": True,
            "warning": None,
        }

    nodes = sorted(G.nodes())
    partitions: list[list[int]] = []

    for i in range(int(n_runs)):
        seed_i = int(base_seed) + i
        try:
            communities = nx.community.louvain_communities(G, weight="weight_abs", seed=seed_i)
        except Exception:
            communities = list(nx.community.greedy_modularity_communities(G, weight="weight_abs"))
        membership = {node: cid for cid, members in enumerate(communities) for node in members}
        partitions.append([membership.get(n, -1) for n in nodes])

    ami_scores: list[float] = []
    for i in range(len(partitions)):
        for j in range(i + 1, len(partitions)):
            ami_scores.append(float(adjusted_mutual_info_score(partitions[i], partitions[j])))

    ami_mean = float(np.mean(ami_scores)) if ami_scores else 1.0
    ami_min = float(np.min(ami_scores)) if ami_scores else 1.0
    stable = ami_mean >= 0.8

    warning = None
    if not stable:
        warning = (
            f"Community assignments show low stability (mean AMI = {ami_mean:.3f} across "
            f"{n_runs} Louvain runs). Consider adjusting CORR_ABS_THRESHOLD to produce a "
            "cleaner community structure before interpreting module memberships."
        )

    return {
        "ami_mean": round(ami_mean, 4),
        "ami_min": round(ami_min, 4),
        "n_runs": int(n_runs),
        "stable": stable,
        "warning": warning,
    }


# ============================================================
# Rewiring Utilities (optional, but handy for run.py)
# ============================================================

def edge_id(a: str, b: str) -> str:
    """Stable undirected edge key."""
    a = str(a)
    b = str(b)
    return f"{a}__{b}" if a <= b else f"{b}__{a}"


def add_edge_ids(edges_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Return a copy of edges_df with an 'edge_id' column."""
    if edges_df is None or edges_df.empty:
        return pd.DataFrame(columns=EDGE_COLS + ["edge_id"])

    df = edges_df.copy()
    for c in ["node_a", "node_b"]:
        if c not in df.columns:
            raise ValueError(f"edges_df missing required column: '{c}'")
    df["node_a"] = df["node_a"].astype(str)
    df["node_b"] = df["node_b"].astype(str)
    df["edge_id"] = [edge_id(a, b) for a, b in zip(df["node_a"], df["node_b"])]
    return df


def compute_rewiring(
    edges_a: Optional[pd.DataFrame],
    edges_b: Optional[pd.DataFrame],
    *,
    label_a: str,
    label_b: str,
    sign_switch_min_r: float = 0.0,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute rewiring tables between two edge sets.

    Parameters
    ----------
    edges_a, edges_b
        Edge DataFrames with columns: node_a, node_b, r, q, sign.
    label_a, label_b
        Group labels used to name per-condition columns.
    sign_switch_min_r
        Minimum |r| required in BOTH networks for an edge to be flagged as a sign switch.
        Setting this > 0 (e.g. 0.5) prevents marginal correlations near zero from being
        counted as sign switches due to sampling noise.  Default 0.0 (no filter).

    Returns
    -------
    shared, only_a, only_b, union_with_stats
        union_with_stats includes per-condition r/q/sign where available and a 'state' column:
        shared / <label_a>_only / <label_b>_only
    """
    ea = add_edge_ids(edges_a)
    eb = add_edge_ids(edges_b)

    A = set(ea["edge_id"].tolist())
    B = set(eb["edge_id"].tolist())

    shared_ids = A & B
    only_a_ids = A - B
    only_b_ids = B - A

    shared = ea[ea["edge_id"].isin(shared_ids)].copy()
    only_a = ea[ea["edge_id"].isin(only_a_ids)].copy()
    only_b = eb[eb["edge_id"].isin(only_b_ids)].copy()

    # Prepare per-condition columns for union
    ea2 = ea.rename(
        columns={"r": f"r_{label_a}", "q": f"q_{label_a}", "sign": f"sign_{label_a}"}
    )
    eb2 = eb.rename(
        columns={"r": f"r_{label_b}", "q": f"q_{label_b}", "sign": f"sign_{label_b}"}
    )

    union = pd.merge(
        ea2[["edge_id", "node_a", "node_b", f"r_{label_a}", f"q_{label_a}", f"sign_{label_a}"]],
        eb2[["edge_id", "node_a", "node_b", f"r_{label_b}", f"q_{label_b}", f"sign_{label_b}"]],
        on="edge_id",
        how="outer",
        suffixes=("", "_dup"),
    )

    # Resolve possible duplicates from merge
    for col in ["node_a", "node_b"]:
        dup = f"{col}_dup"
        if dup in union.columns:
            union[col] = union[col].fillna(union[dup])
            union = union.drop(columns=[dup])

    union[f"present_{label_a}"] = union[f"r_{label_a}"].notna()
    union[f"present_{label_b}"] = union[f"r_{label_b}"].notna()

    def _state(row) -> str:
        pa = bool(row[f"present_{label_a}"])
        pb = bool(row[f"present_{label_b}"])
        if pa and pb:
            return "shared"
        if pa and not pb:
            return f"{label_a}_only"
        if pb and not pa:
            return f"{label_b}_only"
        return "absent"

    union["state"] = union.apply(_state, axis=1)

    _min_r = float(sign_switch_min_r)

    def _sign_switch(row) -> bool:
        if row["state"] != "shared":
            return False
        sa = row.get(f"sign_{label_a}", None)
        sb = row.get(f"sign_{label_b}", None)
        if pd.isna(sa) or pd.isna(sb):
            return False
        if _min_r > 0.0:
            ra = row.get(f"r_{label_a}", None)
            rb = row.get(f"r_{label_b}", None)
            if pd.isna(ra) or pd.isna(rb):
                return False
            if abs(float(ra)) < _min_r or abs(float(rb)) < _min_r:
                return False
        return str(sa) != str(sb)

    union["sign_switch"] = union.apply(_sign_switch, axis=1)

    return shared, only_a, only_b, union


def summarize_rewiring(
    union_df: pd.DataFrame,
    *,
    label_a: str,
    label_b: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if union_df is None or union_df.empty:
        empty_summary = pd.DataFrame(
            [
                {
                    "n_shared_edges": 0,
                    f"n_{label_a}_only_edges": 0,
                    f"n_{label_b}_only_edges": 0,
                    "n_sign_switch_edges": 0,
                }
            ]
        )
        empty_nodes = pd.DataFrame(
            columns=[
                "feature_id",
                "shared_edges",
                f"{label_a}_only_edges",
                f"{label_b}_only_edges",
                "sign_switch_edges",
                "rewiring_score",
            ]
        )
        return empty_summary, empty_nodes

    n_shared = int((union_df["state"] == "shared").sum())
    n_a_only = int((union_df["state"] == f"{label_a}_only").sum())
    n_b_only = int((union_df["state"] == f"{label_b}_only").sum())
    n_sign_switch = int(union_df["sign_switch"].fillna(False).sum())
    total_union = n_shared + n_a_only + n_b_only
    total_state_specific = n_a_only + n_b_only
    summary = pd.DataFrame(
        [
            {
                "n_shared_edges": n_shared,
                f"n_{label_a}_only_edges": n_a_only,
                f"n_{label_b}_only_edges": n_b_only,
                "n_sign_switch_edges": n_sign_switch,
                "total_edges_in_union": total_union,
                "total_state_specific_edges": total_state_specific,
                # fraction of union edges that are state-specific (0 = fully conserved, 1 = fully rewired)
                "rewiring_proportion": round(total_state_specific / max(1, total_union), 4),
                # fraction of shared edges that invert sign
                "sign_switch_proportion": round(n_sign_switch / max(1, n_shared), 4),
            }
        ]
    )

    node_scores: Dict[str, Dict[str, float]] = {}
    for row in union_df.itertuples(index=False):
        state = str(row.state)
        sign_switch = bool(getattr(row, "sign_switch", False))
        for node in [str(row.node_a), str(row.node_b)]:
            node_scores.setdefault(
                node,
                {
                    "feature_id": node,
                    "shared_edges": 0,
                    f"{label_a}_only_edges": 0,
                    f"{label_b}_only_edges": 0,
                    "sign_switch_edges": 0,
                },
            )
            if state == "shared":
                node_scores[node]["shared_edges"] += 1
            elif state == f"{label_a}_only":
                node_scores[node][f"{label_a}_only_edges"] += 1
            elif state == f"{label_b}_only":
                node_scores[node][f"{label_b}_only_edges"] += 1
            if sign_switch:
                node_scores[node]["sign_switch_edges"] += 1

    nodes = pd.DataFrame(node_scores.values())
    # Absolute rewiring score (higher = more rewired edges in total)
    nodes["rewiring_score"] = (
        nodes[f"{label_a}_only_edges"].astype(float)
        + nodes[f"{label_b}_only_edges"].astype(float)
        + nodes["sign_switch_edges"].astype(float)
    )
    # Normalized rewiring fraction: state-specific / total union degree (0=conserved, 1=fully rewired)
    nodes["total_union_degree"] = (
        nodes["shared_edges"].astype(float)
        + nodes[f"{label_a}_only_edges"].astype(float)
        + nodes[f"{label_b}_only_edges"].astype(float)
    )
    nodes["rewiring_fraction"] = (
        (nodes[f"{label_a}_only_edges"].astype(float) + nodes[f"{label_b}_only_edges"].astype(float))
        / nodes["total_union_degree"].clip(lower=1.0)
    ).round(4)
    nodes = nodes.sort_values(
        ["rewiring_score", "sign_switch_edges", f"{label_a}_only_edges", f"{label_b}_only_edges", "feature_id"],
        ascending=[False, False, False, False, True],
    )
    return summary, nodes


# ============================================================
# Export Network
# ============================================================

def export_network(
    selected_features: pd.DataFrame,
    edges_df: Optional[pd.DataFrame],
    outdir: Path,
) -> None:
    """
    Export a Cytoscape-ready network:
      - nodes.csv
      - edges.csv
      - network.graphml

    Notes
    -----
    - Node identifier is 'feature_id' (string), derived from selected_features.index.
    - GraphML attributes are cast to scalar/string to avoid serialization issues.
    - This function does NOT add pathway columns; those can be merged later in Cytoscape.
    """

    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Nodes table
    # --------------------------------------------------------
    nodes_df = selected_features.copy()

    # Unique stable ID = row ID (index)
    nodes_df.insert(0, "feature_id", nodes_df.index.astype(str))

    # Ensure readable name column exists
    if "name" not in nodes_df.columns:
        nodes_df["name"] = nodes_df.index.astype(str)

    # GraphML-safe casting
    for col in nodes_df.columns:
        if col == "feature_id":
            continue
        nodes_df[col] = nodes_df[col].apply(
            lambda x: x
            if np.isscalar(x) and not pd.isna(x)
            else ("" if pd.isna(x) else str(x))
        )

    nodes_df.to_csv(outdir / "nodes.csv", index=False)
    metric_cols = [col for col in NODE_METRIC_COLS if col in nodes_df.columns]
    if metric_cols:
        metric_export_cols = ["feature_id"]
        if "name" in nodes_df.columns:
            metric_export_cols.append("name")
        metric_export_cols.extend(metric_cols)
        nodes_df[metric_export_cols].to_csv(outdir / "node_metrics.csv", index=False)

    # --------------------------------------------------------
    # Edges table
    # --------------------------------------------------------
    if edges_df is None or edges_df.empty:
        edges_out = pd.DataFrame(columns=EDGE_COLS)
    else:
        edges_out = edges_df.copy()
        for c in ["node_a", "node_b", "r", "q"]:
            if c not in edges_out.columns:
                raise ValueError(f"edges_df missing required column: '{c}'")

        edges_out["node_a"] = edges_out["node_a"].astype(str)
        edges_out["node_b"] = edges_out["node_b"].astype(str)
        edges_out["r"] = edges_out["r"].astype(float)
        edges_out["q"] = edges_out["q"].astype(float)

        if "sign" not in edges_out.columns:
            edges_out["sign"] = np.where(edges_out["r"] >= 0, "pos", "neg")
        else:
            edges_out["sign"] = edges_out["sign"].astype(str)

        edges_out = edges_out[EDGE_COLS]

    edges_out.to_csv(outdir / "edges.csv", index=False)

    # --------------------------------------------------------
    # GraphML export
    # --------------------------------------------------------
    G = nx.Graph()

    # Add nodes
    for _, row in nodes_df.iterrows():
        fid = str(row["feature_id"])
        attrs = row.to_dict()
        G.add_node(fid, **attrs)

    # Add edges
    if not edges_out.empty:
        for a, b, r, q, sign in edges_out.itertuples(index=False):
            G.add_edge(str(a), str(b), r=float(r), q=float(q), sign=str(sign))

    nx.write_graphml(G, outdir / "network.graphml")

    logger.info(
        "Network exported to %s with %d nodes and %d edges.",
        str(outdir),
        G.number_of_nodes(),
        G.number_of_edges(),
    )


def export_rewiring_network(
    node_annotations: pd.DataFrame,
    union_df: pd.DataFrame,
    rewiring_nodes: pd.DataFrame,
    outpath: Path,
) -> None:
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    G = nx.Graph()

    nodes_df = node_annotations.copy()
    nodes_df.index = nodes_df.index.astype(str)
    if rewiring_nodes is not None and not rewiring_nodes.empty:
        rew_nodes = rewiring_nodes.copy()
        rew_nodes["feature_id"] = rew_nodes["feature_id"].astype(str)
        rew_nodes = rew_nodes.set_index("feature_id")
        nodes_df = nodes_df.join(rew_nodes, how="left")

    nodes_df["feature_id"] = nodes_df.index.astype(str)
    for _, row in nodes_df.iterrows():
        attrs = {}
        for key, value in row.to_dict().items():
            if pd.isna(value):
                attrs[key] = ""
            elif np.isscalar(value):
                attrs[key] = value
            else:
                attrs[key] = str(value)
        G.add_node(str(row["feature_id"]), **attrs)

    if union_df is not None and not union_df.empty:
        for _, row in union_df.iterrows():
            attrs = {}
            for key, value in row.to_dict().items():
                if key in {"node_a", "node_b", "edge_id"}:
                    continue
                if pd.isna(value):
                    attrs[key] = ""
                elif np.isscalar(value):
                    attrs[key] = value
                else:
                    attrs[key] = str(value)
            G.add_edge(str(row["node_a"]), str(row["node_b"]), **attrs)

    nx.write_graphml(G, outpath)
    logger.info(
        "Rewiring network exported to %s with %d nodes and %d edges.",
        str(outpath),
        G.number_of_nodes(),
        G.number_of_edges(),
    )
