"""
Method comparison benchmark for PhenoRewire.

Compares three approaches for identifying rewired/differentially-correlated nodes
on a shared synthetic ground-truth dataset:

1. **PhenoRewire**      — full pipeline (feature selection + network + rewiring)
2. **Naive baseline**   — top nodes by |log2FC| of mean intensities between groups
3. **DiffCorr-style**   — Fisher Z-test on pairwise correlations (r_ref vs r_case),
                          nodes scored by max -log10(p) across all their edge tests
                          (pure Python reimplementation, no R dependency)

All three are evaluated with:
  - AUROC    over the 8 injected rewired nodes (F000-F007)
  - AP       average precision
  - P@8      precision at 8 (exact ground-truth size)

Synthetic dataset is identical to the rewiring benchmark:
  N_NODES=20, N_REWIRED=8, n_per_group in [10, 20, 30, 50].

Usage::

    python -m benchmarks.comparison.method_comparison   # from repo root
    python benchmarks/comparison/method_comparison.py   # direct

Output:
    - Prints a comparison table (n x method x metric)
    - Saves method_comparison_results.csv next to this script
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score
from statsmodels.stats.multitest import multipletests

from phenorewire.networks import correlation_network, compute_rewiring, summarize_rewiring
from phenorewire.stats import fast_spearman_permutation_test, adaptive_fdr_selection

# ---------------------------------------------------------------------------
# Shared constants (identical to rewiring_benchmark.py)
# ---------------------------------------------------------------------------

N_NODES = 20
N_REWIRED = 8
SEED = 42
SAMPLE_SIZES = [10, 20, 30, 50]


# ---------------------------------------------------------------------------
# Synthetic data (copy of rewiring_benchmark make_synthetic_dataset)
# ---------------------------------------------------------------------------

def _pos_def_cov(n: int, off_diag: float) -> np.ndarray:
    C = np.full((n, n), off_diag)
    np.fill_diagonal(C, 1.0)
    C += np.eye(n) * 0.05
    return C


def make_synthetic_dataset(
    n_per_group: int,
    seed: int = SEED,
) -> tuple[pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(seed)
    feature_ids = [f"F{i:03d}" for i in range(N_NODES)]
    sample_ids = [f"S{i:03d}" for i in range(2 * n_per_group)]

    C_ref = _pos_def_cov(N_NODES, 0.0)
    for i in range(N_REWIRED):
        for j in range(N_REWIRED):
            if i != j:
                C_ref[i, j] = 0.75

    C_case = _pos_def_cov(N_NODES, 0.0)
    for i in range(N_REWIRED):
        for j in range(N_REWIRED):
            if i != j:
                if (i, j) in {(0, 1), (1, 0), (2, 3), (3, 2)}:
                    C_case[i, j] = -0.72
                else:
                    C_case[i, j] = 0.08

    L_ref = np.linalg.cholesky(C_ref)
    L_case = np.linalg.cholesky(C_case)
    Z_ref = rng.standard_normal((n_per_group, N_NODES)) @ L_ref.T
    Z_case = rng.standard_normal((n_per_group, N_NODES)) @ L_case.T
    Z_case[:, :N_REWIRED] += 1.2

    X_all = np.vstack([Z_ref, Z_case]).T
    X = pd.DataFrame(X_all, index=feature_ids, columns=sample_ids)
    y = np.array([0] * n_per_group + [1] * n_per_group, dtype=int)
    return X, y


GT_NODES = {f"F{i:03d}" for i in range(N_REWIRED)}


# ---------------------------------------------------------------------------
# Method 1 — PhenoRewire
# ---------------------------------------------------------------------------

def run_phenorewire(X: pd.DataFrame, y: np.ndarray, seed: int) -> pd.Series:
    """Return a Series indexed by feature_id with rewiring_score (higher = more rewired)."""
    mat = X.values.astype(float)

    r_obs, p_emp = fast_spearman_permutation_test(mat, y, n_permutations=200, seed=seed)
    _, qvals, _, _ = multipletests(p_emp, alpha=0.2, method="fdr_bh")
    rej, _ = adaptive_fdr_selection(qvals, base_alpha=0.2, adaptive=True, min_selected=5, alpha_ceiling=0.4)

    selected_ids = X.index[rej].tolist()
    if len(selected_ids) < 3:
        order = np.argsort(-np.abs(np.where(np.isfinite(r_obs), r_obs, 0.0)))
        selected_ids = X.index[order[:max(5, N_REWIRED)]].tolist()

    samp_ids = X.columns.tolist()
    ref_samps = [samp_ids[i] for i in range(X.shape[1]) if y[i] == 0]
    case_samps = [samp_ids[i] for i in range(X.shape[1]) if y[i] == 1]

    X_sel = X.loc[selected_ids]
    X_ref = X_sel[ref_samps]
    X_case = X_sel[case_samps]

    _, edges_ref, _ = correlation_network(
        X_ref, abs_r_thr=0.4, fdr_alpha=0.2,
        adaptive_thresholds=True, min_edges=3, max_density=0.6,
        floor_abs_r_thr=0.3, ceil_fdr_alpha=0.4,
    )
    _, edges_case, _ = correlation_network(
        X_case, abs_r_thr=0.4, fdr_alpha=0.2,
        adaptive_thresholds=True, min_edges=3, max_density=0.6,
        floor_abs_r_thr=0.3, ceil_fdr_alpha=0.4,
    )

    if edges_ref is None or edges_case is None:
        return pd.Series({fid: 0.0 for fid in X.index})

    _, _, _, union = compute_rewiring(edges_ref, edges_case, label_a="ref", label_b="case", sign_switch_min_r=0.3)
    _, rewiring_nodes = summarize_rewiring(union, label_a="ref", label_b="case")

    if rewiring_nodes is None or rewiring_nodes.empty:
        return pd.Series({fid: 0.0 for fid in X.index})

    scores = rewiring_nodes.set_index("feature_id")["rewiring_score"]
    all_scores = pd.Series(0.0, index=X.index)
    all_scores.update(scores)
    return all_scores


# ---------------------------------------------------------------------------
# Method 2 — Naive log2FC
# ---------------------------------------------------------------------------

def run_naive_log2fc(X: pd.DataFrame, y: np.ndarray) -> pd.Series:
    """Score nodes by |log2FC| of group means (shifted to positive first)."""
    mat = X.values.astype(float)
    # Shift to positive (avoids log of negative)
    mat = mat - mat.min(axis=1, keepdims=True) + 1.0

    ref_idx = np.where(y == 0)[0]
    case_idx = np.where(y == 1)[0]

    mean_ref = mat[:, ref_idx].mean(axis=1)
    mean_case = mat[:, case_idx].mean(axis=1)

    log2fc = np.log2(mean_case + 1e-9) - np.log2(mean_ref + 1e-9)
    return pd.Series(np.abs(log2fc), index=X.index)


# ---------------------------------------------------------------------------
# Method 3 — DiffCorr-style (Fisher Z-test)
# ---------------------------------------------------------------------------

def _fisher_z_test(r1: float, r2: float, n1: int, n2: int) -> float:
    """Two-sample Fisher Z-test for difference in Pearson/Spearman correlations.
    Returns two-tailed p-value."""
    # Clamp to avoid atanh blow-up
    r1 = np.clip(r1, -0.9999, 0.9999)
    r2 = np.clip(r2, -0.9999, 0.9999)
    z1 = np.arctanh(r1)
    z2 = np.arctanh(r2)
    se = np.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    z_stat = (z1 - z2) / se
    p = 2.0 * stats.norm.sf(abs(z_stat))
    return float(p)


def run_diffcorr(X: pd.DataFrame, y: np.ndarray) -> pd.Series:
    """
    DiffCorr-style: for each pair of features compute Fisher Z-test on
    Spearman r_ref vs r_case; score each node by its max -log10(p) across
    all incident edge tests.
    """
    ref_idx = np.where(y == 0)[0]
    case_idx = np.where(y == 1)[0]
    n_ref = len(ref_idx)
    n_case = len(case_idx)

    mat = X.values.astype(float)
    mat_ref = mat[:, ref_idx]
    mat_case = mat[:, case_idx]

    n_feat = mat.shape[0]

    # Compute Spearman r matrices
    def spearman_matrix(m: np.ndarray) -> np.ndarray:
        """Return (n_feat x n_feat) Spearman correlation matrix."""
        # Rank each row
        from scipy.stats import rankdata
        ranked = np.array([rankdata(row) for row in m])
        r_mat = np.corrcoef(ranked)
        return r_mat

    r_ref = spearman_matrix(mat_ref)
    r_case = spearman_matrix(mat_case)

    # Node score = max -log10(p) across all incident edges
    node_scores = np.zeros(n_feat)
    for i in range(n_feat):
        max_neg_log_p = 0.0
        for j in range(i + 1, n_feat):
            p = _fisher_z_test(r_ref[i, j], r_case[i, j], n_ref, n_case)
            neg_log_p = -np.log10(p + 1e-300)
            if neg_log_p > max_neg_log_p:
                max_neg_log_p = neg_log_p
        node_scores[i] = max_neg_log_p

    # Second pass: also update node j from pairs (i,j) where j > i
    node_scores2 = np.zeros(n_feat)
    for i in range(n_feat):
        for j in range(i + 1, n_feat):
            p = _fisher_z_test(r_ref[i, j], r_case[i, j], n_ref, n_case)
            neg_log_p = -np.log10(p + 1e-300)
            if neg_log_p > node_scores2[i]:
                node_scores2[i] = neg_log_p
            if neg_log_p > node_scores2[j]:
                node_scores2[j] = neg_log_p

    return pd.Series(node_scores2, index=X.index)


# ---------------------------------------------------------------------------
# Recovery metrics
# ---------------------------------------------------------------------------

def recovery_metrics(scores: pd.Series, gt_nodes: set[str]) -> dict[str, float]:
    ids = scores.index.tolist()
    sc = scores.fillna(0.0).tolist()
    binary = [1 if fid in gt_nodes else 0 for fid in ids]

    if sum(binary) == 0 or sum(binary) == len(binary):
        return {"auroc": float("nan"), "ap": float("nan"), "p_at_k": float("nan")}

    auroc = float(roc_auc_score(binary, sc))
    ap = float(average_precision_score(binary, sc))

    k = len(gt_nodes)
    top_k = scores.nlargest(k).index.tolist()
    p_at_k = sum(1 for f in top_k if f in gt_nodes) / k

    return {
        "auroc": round(auroc, 4),
        "ap": round(ap, 4),
        "p_at_k": round(p_at_k, 4),
    }


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_comparison(seeds: list[int] | None = None) -> pd.DataFrame:
    if seeds is None:
        seeds = [SEED, SEED + 1, SEED + 2]

    methods = {
        "PhenoRewire": lambda X, y, s: run_phenorewire(X, y, s),
        "Naive_log2FC": lambda X, y, s: run_naive_log2fc(X, y),
        "DiffCorr_Fisher": lambda X, y, s: run_diffcorr(X, y),
    }

    records = []
    for n in SAMPLE_SIZES:
        for method_name, method_fn in methods.items():
            metrics_per_seed: list[dict[str, float]] = []
            for s in seeds:
                X, y = make_synthetic_dataset(n_per_group=n, seed=s)
                scores = method_fn(X, y, s)
                m = recovery_metrics(scores, GT_NODES)
                metrics_per_seed.append(m)

            avg: dict[str, float] = {}
            for k in metrics_per_seed[0]:
                vals = [m[k] for m in metrics_per_seed if not np.isnan(m.get(k, float("nan")))]
                avg[k] = round(float(np.mean(vals)), 4) if vals else float("nan")

            records.append({
                "n_per_group": n,
                "n_total": 2 * n,
                "method": method_name,
                **avg,
            })

    return pd.DataFrame(records)


def main() -> None:
    print("PhenoRewire -- Method Comparison Benchmark")
    print("=" * 56)
    print(f"Features: {N_NODES}  |  Ground-truth rewired nodes: {N_REWIRED} (F000-F007)")
    print(f"Sample sizes: {SAMPLE_SIZES}  |  Seeds averaged: 3")
    print()

    results = run_comparison()

    # Pretty-print as a pivot: rows = n, columns = method x metric
    pivot = results.pivot(index="n_per_group", columns="method", values=["auroc", "ap", "p_at_k"])
    pivot.columns = [f"{m}_{c}" for c, m in pivot.columns]
    pivot = pivot.reset_index()
    print(pivot.to_string(index=False))

    out_path = Path(__file__).parent / "method_comparison_results.csv"
    results.to_csv(out_path, index=False)
    print(f"\nFull results saved to: {out_path}")

    # Summary: best method at each n
    print("\nBest AUROC per sample size:")
    for n, grp in results.groupby("n_per_group"):
        best = grp.loc[grp["auroc"].idxmax()]
        print(f"  n={int(n):3d}/group  ->  {best['method']:20s}  AUROC={best['auroc']:.4f}")


if __name__ == "__main__":
    main()
