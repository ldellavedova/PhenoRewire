# ── PhenoRewire · Statistical utilities ──────────────────────────────────────
# Empirical permutation testing, adaptive FDR threshold relaxation, and
# association statistics. Called by phenotype_selection.py and temporal_selection.py.
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.stats import rankdata, spearmanr

# Tolerance for the empirical p-value comparison.  With a discrete target (e.g. a
# binary phenotype) many permuted correlations are *exactly* equal to the observed
# one; comparing them with a bare ">=" makes the count depend on last-ulp floating
# point noise, which differs across BLAS builds and SciPy versions.  Counting
# near-equal values as ties keeps empirical p-values reproducible.
_TIE_TOL = 1e-12

# Upper bound on the permutation matrix materialised at once, in floats.
_CHUNK_CELLS = 4_000_000


def _rank_normalize(a: np.ndarray) -> np.ndarray:
    """Center and L2-normalize rank rows so that ``a @ b`` is a Spearman rho.

    Rows with zero variance (all values tied) become all-zero rows, which yield
    rho = 0 downstream rather than a division by zero.
    """
    a = np.atleast_2d(np.asarray(a, dtype=float))
    centered = a - a.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    return np.divide(centered, norms, out=np.zeros_like(centered), where=norms > 0)


def _permute_one_feature(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int,
    seed_i: int,
) -> tuple[float, float]:
    """Observed Spearman r and empirical p-value for a single feature row.

    Fallback path, used only for rows containing NaN: pairwise deletion changes
    the sample set per feature, so the vectorized form does not apply.
    """
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

    p = (float(np.sum(r_perm >= abs(float(r)) - _TIE_TOL)) + 1.0) / (float(n_perm) + 1.0)
    return (float(r), p)


def spearman_permutation_test(
    mat: np.ndarray,
    target: np.ndarray,
    *,
    n_permutations: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute observed Spearman correlations and empirical p-values per feature.

    Spearman's rho is Pearson's r on ranks, so ranking once up front turns each
    permutation into a matrix product instead of a fresh ``spearmanr`` call.
    Permuting the target is equivalent to permuting its ranks (a permutation does
    not change the tie structure), so the null distribution is unchanged.

    Each feature draws from its own generator seeded with
    ``seed_i = (seed + feature_index) % 2**31``, so results are reproducible and
    independent of how features are chunked.

    Rows containing NaN fall back to a per-feature ``spearmanr`` loop with
    pairwise deletion.  Constant rows return ``(nan, 1.0)``.

    Returns
    -------
    (r_obs, p_emp)
        Observed rho per feature and its empirical two-sided p-value.
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
    if np.isnan(y).any():
        raise ValueError("target contains NaN; drop or impute those samples first.")

    n_features, n_samples = X.shape
    r_obs = np.full(n_features, np.nan, dtype=float)
    p_emp = np.ones(n_features, dtype=float)

    # Ranks of the target, centered and normalized once.
    y_rank_norm = _rank_normalize(rankdata(y))[0]

    with np.errstate(invalid="ignore"):
        varying_rows = np.nanstd(X, axis=1) > 0
    complete_rows = ~np.isnan(X).any(axis=1)
    fast_rows = np.flatnonzero(complete_rows & varying_rows)
    slow_rows = np.flatnonzero(~complete_rows & varying_rows)

    if fast_rows.size:
        X_rank_norm = _rank_normalize(np.apply_along_axis(rankdata, 1, X[fast_rows]))
        r_fast = X_rank_norm @ y_rank_norm
        r_obs[fast_rows] = r_fast

        # Permutations are drawn per feature (to keep the per-feature seed) but the
        # null correlations are computed as one matrix product per feature.
        block = max(1, min(n_perm, _CHUNK_CELLS // max(n_samples, 1)))
        for local_idx, feature_idx in enumerate(fast_rows):
            rng = np.random.default_rng((int(seed) + int(feature_idx)) % (2 ** 31))
            threshold = abs(float(r_fast[local_idx])) - _TIE_TOL
            x_row = X_rank_norm[local_idx]
            n_ge = 0
            drawn = 0
            while drawn < n_perm:
                size = min(block, n_perm - drawn)
                perms = np.empty((size, n_samples), dtype=float)
                for j in range(size):
                    perms[j] = rng.permutation(y_rank_norm)
                n_ge += int(np.count_nonzero(np.abs(perms @ x_row) >= threshold))
                drawn += size
            p_emp[int(feature_idx)] = (float(n_ge) + 1.0) / (float(n_perm) + 1.0)

    for feature_idx in slow_rows:
        idx = int(feature_idx)
        r_obs[idx], p_emp[idx] = _permute_one_feature(
            X[idx], y, n_perm, (int(seed) + idx) % (2 ** 31)
        )

    p_emp = np.where(np.isfinite(p_emp), p_emp, 1.0)
    return r_obs, p_emp


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
