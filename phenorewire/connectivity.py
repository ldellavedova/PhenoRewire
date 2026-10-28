# ── PhenoRewire · Step 2b — Connectivity-based feature selection ─────────────
# Selects features whose correlation neighbourhood changes between two states,
# whatever their abundance does.  This reaches what abundance-based selection
# cannot: a metabolite that keeps its level but changes its partners.
#
# Score.  For feature i, its Spearman correlations to every other feature are
# computed once per state; the score is the mean of the TOP_K largest absolute
# changes.  The top-k mean, rather than the mean over all partners, keeps the
# score sensitive to a feature that rewires against a handful of partners
# (validated on simulated ground truth: dev/benchmarks/null_benchmark.py).
#
# Significance.  State labels are permuted across the pooled samples and the
# score recomputed.  Unlike the rewiring null (nulls.py), features are NOT
# centred within state first: on simulated ground truth both versions have the
# same power, but the uncentred null kept the empirical FDR at 3.7-7.2% across
# scenarios and sample sizes, against 4.3-7.8% centred, and the difference was
# largest exactly when abundance also shifts (3.7-3.9% vs 7.3-7.8%).  For the
# rewiring null the opposite holds: uncentred, it is invalid.  A per-feature
# empirical p is bounded
# below by 1 / (n_permutations + 1), which BH over hundreds of features cannot
# clear at any practical permutation count, so each feature is standardised
# against its own null and its p-value is read off the pooled distribution of
# standardised null scores.  The per-feature empirical p is reported alongside.
from __future__ import annotations

from typing import Any

import logging

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

# Below this many samples in the smaller state the selector recovers almost
# nothing in simulation (dev/benchmarks), so it is not run.
CONNECTIVITY_MIN_N = 14


def spearman_matrix(mat: np.ndarray) -> np.ndarray:
    """Spearman correlation between the rows of ``mat`` (features x samples).

    Pearson on ranks.  Constant rows give NaN rather than collapsing the whole
    result, as ``scipy.stats.spearmanr`` can.
    """
    ranks = rankdata(np.asarray(mat, dtype=float), axis=1)
    centred = ranks - ranks.mean(axis=1, keepdims=True)
    sd = np.sqrt((centred ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (centred @ centred.T) / np.outer(sd, sd)
    out[~np.isfinite(out)] = np.nan
    return out


def differential_connectivity(a: np.ndarray, b: np.ndarray, top_k: int) -> np.ndarray:
    """Per feature, the mean of its ``top_k`` largest |rho_a - rho_b| over partners."""
    d = np.abs(spearman_matrix(a) - spearman_matrix(b))
    np.fill_diagonal(d, np.nan)
    d = np.nan_to_num(d, nan=-np.inf)
    k = min(int(top_k), max(d.shape[1] - 1, 1))
    # The k largest per row, without sorting the whole row.
    top = np.partition(d, d.shape[1] - k, axis=1)[:, -k:]
    with np.errstate(invalid="ignore"):
        out = top.mean(axis=1)
    out[~np.isfinite(out)] = np.nan
    return out


def connectivity_selection(
    X_a: pd.DataFrame,
    X_b: pd.DataFrame,
    *,
    top_k: int,
    n_permutations: int,
    fdr_alpha: float,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Score every feature for change in its correlation neighbourhood.

    ``X_a`` and ``X_b`` are features x samples for the two states, with the
    same features in the same order.  Returns a per-feature table and a summary.
    """
    feats = [str(f) for f in X_a.index]
    if [str(f) for f in X_b.index] != feats:
        raise ValueError("X_a and X_b must hold the same features in the same order.")
    A = np.asarray(X_a.values, dtype=float)
    B = np.asarray(X_b.values, dtype=float)
    A = np.where(np.isfinite(A), A, 0.0)
    B = np.where(np.isfinite(B), B, 0.0)

    # A feature constant in both states has no correlation neighbourhood at all.
    usable = (np.std(A, axis=1) > 0) | (np.std(B, axis=1) > 0)
    obs = differential_connectivity(A, B, top_k)

    pooled = np.hstack([A, B])
    n_a = A.shape[1]
    rng = np.random.default_rng(int(seed))
    null = np.empty((int(n_permutations), len(feats)))
    for i in range(null.shape[0]):
        idx = rng.permutation(pooled.shape[1])
        null[i] = differential_connectivity(pooled[:, idx[:n_a]], pooled[:, idx[n_a:]], top_k)

    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nanmean(null, axis=0)
        sd = np.nanstd(null, axis=0, ddof=1)
        z_obs = (obs - mu) / sd
        z_null = (null - mu) / sd
    flat = np.sort(z_null[np.isfinite(z_null)])
    p_pooled = np.ones(len(feats))
    if flat.size:
        finite = np.isfinite(z_obs)
        n_ge = flat.size - np.searchsorted(flat, z_obs[finite], side="left")
        p_pooled[finite] = (n_ge + 1.0) / (flat.size + 1.0)
    ok = np.isfinite(null) & np.isfinite(obs)[None, :]
    p_feature = (((null >= obs[None, :]) & ok).sum(axis=0) + 1.0) / (ok.sum(axis=0) + 1.0)
    p_pooled[~usable | ~np.isfinite(obs)] = 1.0
    p_feature[~usable | ~np.isfinite(obs)] = 1.0

    q = multipletests(p_pooled, method="fdr_bh")[1]
    selected = (q <= float(fdr_alpha)) & usable

    table = pd.DataFrame(
        {
            "feature_id": feats,
            "connectivity_score": obs,
            "connectivity_null_mean": mu,
            "connectivity_z": z_obs,
            "p_connectivity": p_pooled,
            "p_connectivity_per_feature": p_feature,
            "q_connectivity": q,
            "selected_connectivity": selected,
        }
    ).set_index("feature_id")
    summary = {
        "n_features_total": len(feats),
        "n_features_selected": int(selected.sum()),
        "n_a": int(A.shape[1]),
        "n_b": int(B.shape[1]),
        "top_k": int(top_k),
        "n_permutations": int(n_permutations),
        "fdr_alpha": float(fdr_alpha),
        "centred_within_state": False,
    }
    logger.info("Connectivity selection: %d of %d features selected (q <= %s).",
                summary["n_features_selected"], len(feats), fdr_alpha)
    return table, summary
