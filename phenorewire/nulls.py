# ── PhenoRewire · Null models ────────────────────────────────────────────────
# Three questions a rewiring result has to answer before it is interpreted:
#
#   1. Does a state network hold more edges than features with no joint
#      structure would produce?                       -> network_edge_null()
#   2. Is its modular organisation more than its degree sequence imposes?
#                                                      -> modularity_null()
#   3. Do the two states differ in correlation structure by more than an
#      arbitrary split of the same samples would?      -> rewiring_label_null()
#
# Question 3 is the one a rewiring claim rests on.  Questions 1 and 2 qualify
# the individual networks.  Every replicate is filtered with the thresholds
# the observed network was built with, through the same code path
# (networks.pairwise_spearman / networks.edge_mask), so observed and null
# networks are directly comparable.
from __future__ import annotations

from typing import Any, Optional

import logging

import networkx as nx
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .networks import edge_mask, pairwise_spearman

logger = logging.getLogger(__name__)


# ============================================================
# Shared helpers
# ============================================================

def edge_index_set(mat: np.ndarray, abs_r_thr: float, fdr_alpha: float) -> set[tuple[int, int]]:
    """Edges passing the network filter, as row-index pairs of ``mat``.

    ``mat`` is features x samples.  Constant or invalid rows are dropped before
    correlating and indices are mapped back to the caller's row order, exactly
    as correlation_network() treats them.
    """
    mat = np.array(mat, dtype=float, copy=True)
    if mat.ndim != 2 or mat.shape[0] < 2 or mat.shape[1] < 3:
        return set()
    mat[~np.isfinite(mat)] = np.nan
    with np.errstate(invalid="ignore"):
        std = np.nanstd(mat, axis=1)
    keep = np.isfinite(std) & (std > 0)
    if int(keep.sum()) < 2:
        return set()
    back = np.flatnonzero(keep)
    stats = pairwise_spearman(mat[keep])
    if stats is None:
        return set()
    iu, rvals, _, qvals = stats
    sel = edge_mask(rvals, qvals, abs_r_thr, fdr_alpha)
    return {(int(back[a]), int(back[b])) for a, b in zip(iu[0][sel], iu[1][sel])}


def _empirical_p(null: np.ndarray, observed: float) -> float:
    """One-sided permutation p-value, (#{null >= observed} + 1) / (n + 1)."""
    return float((int(np.sum(null >= observed)) + 1) / (null.size + 1))


def _ratio(observed: float, null_mean: float) -> Optional[float]:
    return round(float(observed) / null_mean, 4) if null_mean > 0 else None


def _null_summary(null: np.ndarray, observed: int) -> dict[str, Any]:
    null_mean = float(null.mean()) if null.size else 0.0
    return {
        "null_mean": round(null_mean, 3),
        "null_sd": round(float(null.std(ddof=1)), 3) if null.size > 1 else None,
        "null_p95": float(np.percentile(null, 95)) if null.size else None,
        "obs_over_null": _ratio(observed, null_mean),
        "p_empirical": _empirical_p(null, observed) if null.size else None,
    }


# ============================================================
# 1. Edges against feature independence
# ============================================================

def network_edge_null(
    X_state: pd.DataFrame,
    *,
    abs_r_thr: float,
    fdr_alpha: float,
    n_permutations: int,
    seed: int,
) -> dict[str, Any]:
    """Observed edge count of one state network against a feature-independence null.

    Each feature's values are permuted independently across the samples of the
    state.  That destroys every between-feature relationship while keeping each
    feature's own distribution, including its ties and zeros.  A network whose
    edge count sits near this null carries no more structure than chance at the
    applied thresholds.
    """
    mat = np.asarray(X_state.values, dtype=float)
    observed = len(edge_index_set(mat, abs_r_thr, fdr_alpha))

    rng = np.random.default_rng(int(seed))
    null = np.empty(int(n_permutations), dtype=float)
    for i in range(null.size):
        perm = np.empty_like(mat)
        for j in range(mat.shape[0]):
            perm[j] = rng.permutation(mat[j])
        null[i] = len(edge_index_set(perm, abs_r_thr, fdr_alpha))

    out = {
        "n_samples": int(mat.shape[1]),
        "n_features": int(mat.shape[0]),
        "applied_abs_r_thr": float(abs_r_thr),
        "applied_fdr_alpha": float(fdr_alpha),
        "observed_edges": int(observed),
    }
    out.update(_null_summary(null, observed))
    out["n_permutations"] = int(null.size)
    return out


# ============================================================
# 2. Modularity against a degree-preserving null
# ============================================================

def _unweighted_graph(edges_df: Optional[pd.DataFrame]) -> nx.Graph:
    G = nx.Graph()
    if edges_df is None or edges_df.empty:
        return G
    G.add_edges_from(zip(edges_df["node_a"].astype(str), edges_df["node_b"].astype(str)))
    return G


def _louvain_q(G: nx.Graph, seed: int) -> float:
    if G.number_of_edges() == 0:
        return 0.0
    communities = nx.community.louvain_communities(G, seed=int(seed))
    return float(nx.community.modularity(G, communities))


def modularity_null(
    edges_df: Optional[pd.DataFrame],
    *,
    n_reps: int,
    seed: int,
) -> dict[str, Any]:
    """Louvain modularity of a network against degree-preserving rewirings of it.

    Raw modularity largely restates density: sparse graphs score high whether or
    not their sparsity means anything.  Comparing against configuration-model
    rewirings, which keep every node's degree, isolates the organisation that
    the degree sequence alone does not impose.  Computed on the unweighted graph
    for both observed and null, so the comparison is like for like; the weighted
    modularity reported in network_summary.csv is unaffected.
    """
    G = _unweighted_graph(edges_df)
    q_obs = _louvain_q(G, seed)
    out: dict[str, Any] = {
        "modularity_unweighted": round(q_obs, 4),
        "modularity_configmodel_mean": None,
        "modularity_configmodel_sd": None,
        "delta_modularity": None,
        "n_configuration_models": 0,
    }
    n_edges = G.number_of_edges()
    if n_edges < 2 or int(n_reps) < 1:
        return out

    rng = np.random.default_rng(int(seed))
    qs = []
    for i in range(int(n_reps)):
        H = G.copy()
        try:
            nx.double_edge_swap(
                H,
                nswap=max(10 * n_edges, 10),
                max_tries=200 * n_edges + 100,
                seed=int(rng.integers(2 ** 31)),
            )
        except (nx.NetworkXError, nx.NetworkXAlgorithmError):
            # Graphs too small or too dense to swap fully keep the swaps made so far.
            pass
        qs.append(_louvain_q(H, int(seed) + i))

    q_cm = float(np.mean(qs))
    out.update(
        {
            "modularity_configmodel_mean": round(q_cm, 4),
            "modularity_configmodel_sd": round(float(np.std(qs, ddof=1)), 4) if len(qs) > 1 else None,
            "delta_modularity": round(q_obs - q_cm, 4),
            "n_configuration_models": len(qs),
        }
    )
    return out


# ============================================================
# 3. State-specific edges against an arbitrary split
# ============================================================

def rewiring_label_null(
    X_a: pd.DataFrame,
    X_b: pd.DataFrame,
    *,
    thresholds_a: tuple[float, float],
    thresholds_b: tuple[float, float],
    n_permutations: int,
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    """State-specific edges of a comparison against a label-permutation null.

    The state labels are permuted across the pooled samples, both networks are
    rebuilt with the thresholds their observed counterparts used, and the edges
    present in exactly one of them are counted.  The observed count against
    that distribution says whether the two states differ in correlation
    structure by more than an arbitrary split of the same samples would.

    Each feature is centred within its own state before pooling.  Without that
    the null is invalid: when the states also differ in abundance, every mixed
    pseudo-group contains samples from both, and two features that both track
    the state correlate strongly inside it.  The null then fills with
    correlations manufactured by the abundance difference.  Centring removes the
    between-state shift and leaves only the question of correlation structure.
    Spearman correlation within a state is unchanged by centring, so the
    observed networks are not affected.

    Returns
    -------
    summary, node_table
        ``summary`` holds the comparison-level result.  ``node_table`` gives,
        for every feature, its observed number of state-specific edges, the null
        mean, an empirical p-value and a BH q-value across features.
    """
    feats = [str(f) for f in X_a.index]
    if [str(f) for f in X_b.index] != feats:
        raise ValueError("X_a and X_b must hold the same features in the same order.")

    A = np.asarray(X_a.values, dtype=float)
    B = np.asarray(X_b.values, dtype=float)
    (r_a, q_a), (r_b, q_b) = thresholds_a, thresholds_b
    n_feat = len(feats)

    def _state_specific(ma: np.ndarray, mb: np.ndarray) -> tuple[int, np.ndarray]:
        diff = edge_index_set(ma, r_a, q_a) ^ edge_index_set(mb, r_b, q_b)
        degree = np.zeros(n_feat, dtype=float)
        for i, j in diff:
            degree[i] += 1
            degree[j] += 1
        return len(diff), degree

    observed, observed_degree = _state_specific(A, B)

    with np.errstate(invalid="ignore"):
        Ac = A - np.nanmean(A, axis=1, keepdims=True)
        Bc = B - np.nanmean(B, axis=1, keepdims=True)
    pooled = np.hstack([Ac, Bc])
    n_a = A.shape[1]

    rng = np.random.default_rng(int(seed))
    null = np.empty(int(n_permutations), dtype=float)
    null_degree = np.empty((null.size, n_feat), dtype=float)
    for i in range(null.size):
        idx = rng.permutation(pooled.shape[1])
        null[i], null_degree[i] = _state_specific(pooled[:, idx[:n_a]], pooled[:, idx[n_a:]])

    summary = {
        "n_a": int(A.shape[1]),
        "n_b": int(B.shape[1]),
        "n_features": n_feat,
        "applied_abs_r_thr_a": float(r_a),
        "applied_fdr_alpha_a": float(q_a),
        "applied_abs_r_thr_b": float(r_b),
        "applied_fdr_alpha_b": float(q_b),
        "observed_state_specific_edges": int(observed),
    }
    summary.update(_null_summary(null, observed))
    summary["n_permutations"] = int(null.size)
    summary["centred_within_state"] = True

    if null.size:
        p_node = (np.sum(null_degree >= observed_degree[None, :], axis=0) + 1.0) / (null.size + 1.0)
        q_node = multipletests(p_node, method="fdr_bh")[1]
        null_mean_node = null_degree.mean(axis=0)
    else:
        p_node = np.full(n_feat, np.nan)
        q_node = np.full(n_feat, np.nan)
        null_mean_node = np.full(n_feat, np.nan)

    node_table = pd.DataFrame(
        {
            "feature_id": feats,
            "observed_state_specific_edges": observed_degree.astype(int),
            "null_mean_state_specific_edges": np.round(null_mean_node, 3),
            "p_empirical": p_node,
            "q_bh": q_node,
        }
    ).sort_values(["p_empirical", "observed_state_specific_edges", "feature_id"], ascending=[True, False, True])
    return summary, node_table.reset_index(drop=True)
