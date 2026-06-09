"""
Ground-truth rewiring benchmark for PhenoRewire.

Builds a synthetic two-group dataset where the rewired edges and sign-switch
edges are injected by design.  Runs the network + rewiring pipeline and
measures how well it recovers the known rewired nodes.

Usage::

    python -m tests.benchmark.rewiring_benchmark      # from repo root
    python tests/benchmark/rewiring_benchmark.py      # direct

Output:
    - Prints a summary table of recovery metrics per sample size
    - Saves rewiring_benchmark_results.csv next to this script

Benchmark definition
--------------------
* N_NODES = 20 features (F000–F019)
* REF group: a random-scale-free-like block correlation structure
* CASE group: 8 "rewired" features (F000–F007) have their edges shuffled:
    - 5 edges that exist in REF are absent in CASE   (REF-only)
    - 5 edges that are absent in REF appear in CASE  (CASE-only)
    - 2 shared edges flip sign                       (sign-switch)
* "Ground-truth rewired nodes": F000–F007 (involved in the injected changes)
* Recovery metric: AUROC of rewiring_score ranking over ground-truth label
* Precision@10: fraction of top-10 ranked nodes that are in ground truth
* Minimum recommended n: first n where AUROC > 0.80
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

from phenorewire.networks import (
    correlation_network,
    compute_rewiring,
    summarize_rewiring,
)
from statsmodels.stats.multitest import multipletests
from phenorewire.stats import fast_spearman_permutation_test, adaptive_fdr_selection

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_NODES = 20
N_REWIRED = 8          # features whose connections change between groups
SEED = 42
SAMPLE_SIZES = [10, 20, 30, 50]


# ---------------------------------------------------------------------------
# Synthetic data generation
# ---------------------------------------------------------------------------

def _pos_def_cov(n: int, off_diag: float) -> np.ndarray:
    C = np.full((n, n), off_diag)
    np.fill_diagonal(C, 1.0)
    C += np.eye(n) * 0.05
    return C


def _make_group_data(
    n_samples: int,
    cov: np.ndarray,
    mean_shift: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return (n_features, n_samples) matrix from multivariate normal."""
    n = cov.shape[0]
    L = np.linalg.cholesky(cov)
    Z = rng.standard_normal((n_samples, n))
    X = Z @ L.T
    X[:, :N_REWIRED] += mean_shift
    return X.T


def make_synthetic_dataset(
    n_per_group: int,
    seed: int = SEED,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, set[str]]]:
    """
    Returns
    -------
    X : pd.DataFrame (N_NODES × 2*n_per_group)
    y : np.ndarray int (0=REF, 1=CASE)
    ground_truth : dict with keys:
        "rewired_nodes"     → set of feature IDs with injected rewiring
        "sign_switch_nodes" → subset of rewired_nodes with injected sign switches
    """
    rng = np.random.default_rng(seed)
    feature_ids = [f"F{i:03d}" for i in range(N_NODES)]
    sample_ids = [f"S{i:03d}" for i in range(2 * n_per_group)]

    # REF covariance: dense block for rewired features
    C_ref = _pos_def_cov(N_NODES, off_diag=0.0)
    # Strong correlations within the rewired block in REF
    for i in range(N_REWIRED):
        for j in range(N_REWIRED):
            if i != j:
                C_ref[i, j] = 0.75

    # CASE covariance: rewired block has different (weaker / sign-flipped) structure
    C_case = _pos_def_cov(N_NODES, off_diag=0.0)
    for i in range(N_REWIRED):
        for j in range(N_REWIRED):
            if i != j:
                # Pairs (0,1) and (2,3) get negative correlations → sign switch
                if (i, j) in {(0, 1), (1, 0), (2, 3), (3, 2)}:
                    C_case[i, j] = -0.72
                else:
                    # Other rewired pairs: near-zero in CASE (lost edge)
                    C_case[i, j] = 0.08

    X_ref = _make_group_data(n_per_group, C_ref, mean_shift=1.2, rng=rng)
    X_case = _make_group_data(n_per_group, C_case, mean_shift=0.0, rng=rng)

    X_all = np.hstack([X_ref, X_case])  # (N_NODES, 2*n_per_group)
    X = pd.DataFrame(X_all, index=feature_ids, columns=sample_ids)
    y = np.array([0] * n_per_group + [1] * n_per_group, dtype=int)

    rewired_nodes = {f"F{i:03d}" for i in range(N_REWIRED)}
    sign_switch_nodes = {"F000", "F001", "F002", "F003"}  # pairs (0,1) and (2,3)

    return X, y, {"rewired_nodes": rewired_nodes, "sign_switch_nodes": sign_switch_nodes}


# ---------------------------------------------------------------------------
# Pipeline (network + rewiring, no file I/O needed)
# ---------------------------------------------------------------------------

def run_pipeline(
    X: pd.DataFrame,
    y: np.ndarray,
    *,
    seed: int = SEED,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """
    Feature selection → REF/CASE networks → rewiring.
    Returns (rewiring_nodes_df, union_df) or (None, None) on failure.
    """
    mat = X.values.astype(float)

    # Feature selection
    r_obs, p_emp = fast_spearman_permutation_test(mat, y, n_permutations=200, seed=seed)
    _, qvals, _, _ = multipletests(p_emp, alpha=0.2, method="fdr_bh")
    rej, _ = adaptive_fdr_selection(
        qvals, base_alpha=0.2, adaptive=True, min_selected=5, alpha_ceiling=0.4
    )

    selected_ids = X.index[rej].tolist()
    if len(selected_ids) < 3:
        order = np.argsort(-np.abs(np.where(np.isfinite(r_obs), r_obs, 0.0)))
        selected_ids = X.index[order[:max(5, N_REWIRED)]].tolist()

    n_samples = X.shape[1]
    n_ref = int((y == 0).sum())
    n_case = int((y == 1).sum())
    samp_ids = X.columns.tolist()
    ref_samps = [samp_ids[i] for i in range(n_samples) if y[i] == 0]
    case_samps = [samp_ids[i] for i in range(n_samples) if y[i] == 1]

    X_sel = X.loc[selected_ids]
    X_ref = X_sel[ref_samps]
    X_case = X_sel[case_samps]

    # Networks
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
        return None, None

    _, _, _, union = compute_rewiring(
        edges_ref, edges_case, label_a="ref", label_b="case", sign_switch_min_r=0.3
    )
    _, rewiring_nodes = summarize_rewiring(union, label_a="ref", label_b="case")

    return rewiring_nodes, union


# ---------------------------------------------------------------------------
# Recovery metrics
# ---------------------------------------------------------------------------

def compute_recovery(
    rewiring_nodes: pd.DataFrame,
    ground_truth: dict[str, set[str]],
    union_df: pd.DataFrame | None,
    all_features: list[str] | None = None,
) -> dict[str, float]:
    """
    Compute recovery metrics.

    Parameters
    ----------
    rewiring_nodes
        Rewiring node summary from summarize_rewiring().
    ground_truth
        Dict with 'rewired_nodes' and 'sign_switch_nodes'.
    union_df
        Union edge DataFrame (used for sign-switch recall).
    all_features
        Full list of feature IDs.  If provided, features absent from
        rewiring_nodes are added with rewiring_score=0 so AUROC has both
        positive and negative examples.
    """
    gt_nodes = ground_truth["rewired_nodes"]
    gt_switch = ground_truth["sign_switch_nodes"]

    if rewiring_nodes is None or rewiring_nodes.empty:
        return {
            "auroc_rewired": float("nan"),
            "ap_rewired": float("nan"),
            "precision_at_10": 0.0,
            "sign_switch_recall": 0.0,
        }

    rew = rewiring_nodes.copy()

    # Add missing features with score=0 so AUROC has negative examples
    if all_features:
        present = set(rew["feature_id"].astype(str).tolist())
        missing = [f for f in all_features if f not in present]
        if missing:
            missing_df = pd.DataFrame({
                "feature_id": missing,
                "rewiring_score": [0.0] * len(missing),
            })
            rew = pd.concat([rew, missing_df], ignore_index=True)

    ids = rew["feature_id"].astype(str).tolist()
    scores = rew["rewiring_score"].fillna(0.0).tolist()

    binary = [1 if fid in gt_nodes else 0 for fid in ids]
    if sum(binary) == 0 or sum(binary) == len(binary):
        auroc = float("nan")
        ap = float("nan")
    else:
        auroc = float(roc_auc_score(binary, scores))
        ap = float(average_precision_score(binary, scores))

    top10 = rewiring_nodes.head(10)["feature_id"].tolist()
    prec10 = sum(1 for f in top10 if f in gt_nodes) / max(len(top10), 1)

    # Sign-switch recall: fraction of known sign-switch nodes that have sign_switch_edges > 0
    switch_recall = 0.0
    if union_df is not None and not union_df.empty and "sign_switch" in union_df.columns:
        switch_edges = union_df[union_df["sign_switch"].fillna(False)]
        detected_switch_nodes = set(switch_edges["node_a"].tolist() + switch_edges["node_b"].tolist())
        if gt_switch:
            switch_recall = len(gt_switch & detected_switch_nodes) / len(gt_switch)

    return {
        "auroc_rewired": round(auroc, 4) if not np.isnan(auroc) else float("nan"),
        "ap_rewired": round(ap, 4) if not np.isnan(ap) else float("nan"),
        "precision_at_10": round(prec10, 4),
        "sign_switch_recall": round(switch_recall, 4),
    }


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_benchmark(seeds: list[int] | None = None) -> pd.DataFrame:
    """Run the benchmark across all sample sizes.  Averages across multiple seeds."""
    if seeds is None:
        seeds = [SEED, SEED + 1, SEED + 2]

    records = []
    for n in SAMPLE_SIZES:
        metrics_per_seed: list[dict[str, float]] = []
        for s in seeds:
            X, y, gt = make_synthetic_dataset(n_per_group=n, seed=s)
            rewiring_nodes, union_df = run_pipeline(X, y, seed=s)
            all_feats = X.index.tolist()
            m = compute_recovery(rewiring_nodes, gt, union_df, all_features=all_feats)
            metrics_per_seed.append(m)

        # Average across seeds
        avg: dict[str, float] = {}
        for k in metrics_per_seed[0]:
            vals = [m[k] for m in metrics_per_seed if not np.isnan(m.get(k, float("nan")))]
            avg[k] = round(float(np.mean(vals)), 4) if vals else float("nan")

        records.append({
            "n_per_group": n,
            "n_total": 2 * n,
            **avg,
        })

    return pd.DataFrame(records)


def main() -> None:
    print("PhenoRewire — Ground-Truth Rewiring Benchmark")
    print("=" * 56)
    print(f"Features: {N_NODES}  |  Rewired nodes (ground truth): {N_REWIRED}")
    print(f"Sample sizes tested: {SAMPLE_SIZES}")
    print()

    results = run_benchmark()
    print(results.to_string(index=False))

    out_path = Path(__file__).parent / "rewiring_benchmark_results.csv"
    results.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")

    # Report minimum reliable n
    reliable = results[results["auroc_rewired"] >= 0.80]
    if not reliable.empty:
        min_n = int(reliable["n_per_group"].min())
        print(f"\nMinimum reliable n per group (AUROC >= 0.80): {min_n}")
        print(f"  -> Recommend >= {min_n} samples per group for rewiring analysis.")
    else:
        print("\nAUROC did not exceed 0.80 at any tested sample size.")
        print("  -> Consider increasing sample size beyond 50 per group.")

    print("\nSign-switch recovery rate:")
    for _, row in results.iterrows():
        print(f"  n={int(row['n_per_group']):3d}/group -> sign_switch_recall = {row['sign_switch_recall']:.3f}")


if __name__ == "__main__":
    main()
