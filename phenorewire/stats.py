# ── PhenoRewire · Statistical utilities ──────────────────────────────────────
# Empirical permutation testing, adaptive FDR threshold relaxation, and
# association statistics. Called by phenotype_selection.py and temporal_selection.py.
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.stats import spearmanr

try:
    from joblib import Parallel, delayed
    _HAS_JOBLIB = True
except ImportError:
    _HAS_JOBLIB = False


def _permute_one_feature(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int,
    seed_i: int,
) -> tuple[float, float]:
    """Compute observed Spearman r and empirical p-value for a single feature row."""
    if np.all(np.isnan(x)) or np.nanstd(x) == 0.0:
        return (np.nan, 1.0)

    r, _ = spearmanr(x, y, nan_policy="omit")
    if not np.isfinite(r):
        return (np.nan, 1.0)

    rng = np.random.default_rng(seed_i)
    r_perm = np.empty(n_perm, dtype=float)
    for j in range(n_perm):
        yp = rng.permutation(y)
        rp, _ = spearmanr(x, yp, nan_policy="omit")
        r_perm[j] = abs(float(rp)) if np.isfinite(rp) else 0.0

    p = (float(np.sum(r_perm >= abs(float(r)))) + 1.0) / (float(n_perm) + 1.0)
    return (float(r), p)


def spearman_permutation_test(
    mat: np.ndarray,
    target: np.ndarray,
    *,
    n_permutations: int,
    seed: int,
    n_jobs: int = -1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute observed Spearman correlations and empirical p-values feature-by-feature.

    Parallelised across features with joblib when available.  Each feature uses a
    deterministic per-feature seed: seed_i = (seed + feature_index) % 2**31, so
    results are fully reproducible regardless of n_jobs.
    """
    X = np.asarray(mat, dtype=float)
    y = np.asarray(target, dtype=float)
    if X.ndim != 2:
        raise ValueError("mat must be 2D (features x samples).")
    if y.ndim != 1 or y.shape[0] != X.shape[1]:
        raise ValueError("target must be 1D and aligned to mat columns.")
    n_perm = int(n_permutations)
    if n_perm < 10:
        raise ValueError("n_permutations must be >= 10 for stability.")

    n_features = X.shape[0]
    global_seed = int(seed)

    seeds = [(global_seed + idx) % (2 ** 31) for idx in range(n_features)]

    if _HAS_JOBLIB and n_jobs != 1:
        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_permute_one_feature)(X[idx], y, n_perm, seeds[idx])
            for idx in range(n_features)
        )
    else:
        results = [
            _permute_one_feature(X[idx], y, n_perm, seeds[idx])
            for idx in range(n_features)
        ]

    r_obs = np.array([r for r, _ in results], dtype=float)
    p_emp = np.array([p for _, p in results], dtype=float)
    p_emp = np.where(np.isfinite(p_emp), p_emp, 1.0)
    return r_obs, p_emp


def fast_spearman_permutation_test(
    mat: np.ndarray,
    target: np.ndarray,
    *,
    n_permutations: int,
    seed: int,
    n_jobs: int = -1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Backward-compatible alias for spearman_permutation_test."""
    return spearman_permutation_test(
        mat,
        target,
        n_permutations=n_permutations,
        seed=seed,
        n_jobs=n_jobs,
    )


def adaptive_fdr_selection(
    qvals: np.ndarray,
    *,
    base_alpha: float,
    adaptive: bool,
    min_selected: int,
    alpha_ceiling: float,
) -> tuple[np.ndarray, dict[str, float | int | bool]]:
    q = np.asarray(qvals, dtype=float)
    finite = np.isfinite(q)
    base_alpha = float(base_alpha)
    alpha_ceiling = float(alpha_ceiling)
    min_selected = int(min_selected)

    base_selected = int((finite & (q <= base_alpha)).sum())
    meta = {
        "base_fdr_alpha": base_alpha,
        "applied_fdr_alpha": base_alpha,
        "adaptive_fdr_used": False,
        "min_selected_target": min_selected,
        "selected_at_base_alpha": base_selected,
        "selected_feature_count": base_selected,
    }
    if not adaptive or base_selected >= min_selected:
        return finite & (q <= base_alpha), meta

    candidate_alphas = sorted(
        {
            round(base_alpha, 6),
            round(alpha_ceiling, 6),
            0.01,
            0.02,
            0.05,
            0.1,
            0.2,
        }
    )
    candidate_alphas = [a for a in candidate_alphas if base_alpha <= a <= alpha_ceiling]
    if not candidate_alphas:
        return finite & (q <= base_alpha), meta

    counts = []
    for alpha in candidate_alphas:
        counts.append((float(alpha), int((finite & (q <= float(alpha))).sum())))

    chosen_alpha = base_alpha
    chosen_count = base_selected
    for alpha, count in counts:
        if count >= min_selected:
            chosen_alpha = alpha
            chosen_count = count
            break
    else:
        chosen_alpha, chosen_count = max(counts, key=lambda item: (item[1], -item[0]))

    meta = {
        "base_fdr_alpha": base_alpha,
        "applied_fdr_alpha": float(chosen_alpha),
        "adaptive_fdr_used": bool(abs(float(chosen_alpha) - base_alpha) > 1e-12),
        "min_selected_target": min_selected,
        "selected_at_base_alpha": base_selected,
        "selected_feature_count": int(chosen_count),
    }
    return finite & (q <= float(chosen_alpha)), meta
