"""
Component-level performance profiling for PhenoRewire.

Times each pipeline stage independently on synthetic matrices of increasing size,
avoiding the statistical constraints (BH-FDR minimum detectable p-value) that make
full end-to-end runs impractical for large synthetic feature sets.

Stages profiled:
  1. Permutation Spearman test     (feature selection)
  2. Correlation network build     (single group)
  3. Two-group rewiring analysis   (compute_rewiring on pre-built networks)

Dataset sizes:
  - 200 features  x  30 samples/group
  - 500 features  x  50 samples/group
  - 1000 features x 100 samples/group  (primary target)

Usage::

    python -m benchmarks.profiling.performance_profile   # from repo root
    python benchmarks/profiling/performance_profile.py   # direct

Output:
    - Prints a timing table per stage and dataset size
    - Saves performance_profile_results.csv next to this script

Notes:
    - Memory tracking requires psutil: pip install psutil
    - The permutation test scales as O(n_features * n_samples * n_permutations)
    - Network construction scales as O(n_selected^2 * n_samples) for pairwise Spearman
    - For n_features >> 200, pre-filter to the top ~50 selected features before
      passing to correlation_network (the pipeline does this automatically)
"""

from __future__ import annotations

import sys
import time
import gc
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

from phenorewire.stats import fast_spearman_permutation_test
from phenorewire.networks import correlation_network, compute_rewiring

# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------

PROFILES = [
    {"n_features": 200,  "n_per_group": 30},
    {"n_features": 500,  "n_per_group": 50},
    {"n_features": 1000, "n_per_group": 100},
]

N_PERMUTATIONS = 200
N_SELECTED = 30   # simulate downstream feature selection output (typical in practice)
SEED = 42


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gen_matrix(n_features: int, n_per_group: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (features x samples) matrix and binary group label y."""
    rng = np.random.default_rng(seed)
    n_signal = min(30, n_features // 5)
    n_total = 2 * n_per_group

    C = np.eye(n_features)
    for i in range(n_signal):
        for j in range(n_signal):
            if i != j:
                C[i, j] = 0.7
    C += np.eye(n_features) * 0.05
    L = np.linalg.cholesky(C)

    Z_ref  = rng.standard_normal((n_per_group, n_features)) @ L.T
    Z_case = rng.standard_normal((n_per_group, n_features)) @ L.T
    Z_case[:, :n_signal] += 1.5

    mat = np.vstack([Z_ref, Z_case]).T  # (n_features, n_total)
    y = np.array([0] * n_per_group + [1] * n_per_group, dtype=int)
    return mat, y


def _gen_selected_df(mat: np.ndarray, feature_ids: list[str], n: int) -> pd.DataFrame:
    """Take the first n rows of mat as a DataFrame (simulates post-selection subset)."""
    sample_ids = [f"S{i:04d}" for i in range(mat.shape[1])]
    return pd.DataFrame(mat[:n, :], index=feature_ids[:n], columns=sample_ids)


def _get_mem_mb() -> float:
    if not HAS_PSUTIL:
        return float("nan")
    return psutil.Process().memory_info().rss / 1024 / 1024


def _time_call(fn, *args, **kwargs) -> tuple[float, object]:
    """Return (elapsed_seconds, result)."""
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    return time.perf_counter() - t0, result


# ---------------------------------------------------------------------------
# Stage timings
# ---------------------------------------------------------------------------

def profile_permutation_test(mat: np.ndarray, y: np.ndarray) -> dict:
    gc.collect()
    t, _ = _time_call(
        fast_spearman_permutation_test, mat, y,
        n_permutations=N_PERMUTATIONS, seed=SEED
    )
    return {"stage": "permutation_test", "t_s": round(t, 3)}


def profile_network_build(mat: np.ndarray, feature_ids: list[str]) -> dict:
    df = _gen_selected_df(mat, feature_ids, N_SELECTED)
    gc.collect()
    t, (_, edges, _) = _time_call(
        correlation_network, df,
        abs_r_thr=0.4, fdr_alpha=0.2,
        adaptive_thresholds=True, min_edges=3, max_density=0.6,
        floor_abs_r_thr=0.3, ceil_fdr_alpha=0.4,
    )
    n_edges = 0 if edges is None else len(edges)
    return {"stage": "network_build", "t_s": round(t, 3), "n_edges": n_edges}


def profile_rewiring(mat: np.ndarray, feature_ids: list[str]) -> dict:
    """Build two slightly different networks and time compute_rewiring."""
    # REF: first half of columns
    df_ref  = _gen_selected_df(mat[:, :mat.shape[1]//2], feature_ids, N_SELECTED)
    df_case = _gen_selected_df(mat[:, mat.shape[1]//2:], feature_ids, N_SELECTED)

    # Build networks (not timed here — part of network_build stage)
    _, edges_ref,  _ = correlation_network(df_ref,  abs_r_thr=0.3, fdr_alpha=0.3,
                                            adaptive_thresholds=True, min_edges=3,
                                            max_density=0.8, floor_abs_r_thr=0.2,
                                            ceil_fdr_alpha=0.5)
    _, edges_case, _ = correlation_network(df_case, abs_r_thr=0.3, fdr_alpha=0.3,
                                            adaptive_thresholds=True, min_edges=3,
                                            max_density=0.8, floor_abs_r_thr=0.2,
                                            ceil_fdr_alpha=0.5)

    if edges_ref is None or edges_case is None:
        return {"stage": "rewiring", "t_s": float("nan"), "note": "empty network"}

    gc.collect()
    t, _ = _time_call(
        compute_rewiring, edges_ref, edges_case,
        label_a="ref", label_b="case", sign_switch_min_r=0.3
    )
    return {"stage": "rewiring", "t_s": round(t, 3)}


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_profiles(seed: int = SEED) -> pd.DataFrame:
    records = []
    for prof in PROFILES:
        n_f = prof["n_features"]
        n_s = prof["n_per_group"]
        label = f"{n_f}f x {n_s}s/g"
        feature_ids = [f"F{i:04d}" for i in range(n_f)]

        print(f"\n  {label}")
        mat, y = _gen_matrix(n_f, n_s, seed)

        stages: list[dict] = []

        print(f"    permutation test ({n_f} features, {N_PERMUTATIONS} perm) ...", end=" ", flush=True)
        r = profile_permutation_test(mat, y)
        print(f"{r['t_s']:.2f}s")
        stages.append(r)

        print(f"    network build ({N_SELECTED} selected features, {2*n_s} samples) ...", end=" ", flush=True)
        r = profile_network_build(mat, feature_ids)
        print(f"{r['t_s']:.2f}s  ({r.get('n_edges', 0)} edges)")
        stages.append(r)

        print(f"    rewiring analysis ({N_SELECTED} features, two groups) ...", end=" ", flush=True)
        r = profile_rewiring(mat, feature_ids)
        print(f"{r.get('t_s', 'N/A')}s")
        stages.append(r)

        total = sum(s.get("t_s", 0.0) or 0.0 for s in stages)
        for s in stages:
            records.append({
                "dataset": label,
                "n_features": n_f,
                "n_per_group": n_s,
                "n_total": 2 * n_s,
                "stage": s["stage"],
                "t_s": s.get("t_s", float("nan")),
            })

        records.append({
            "dataset": label,
            "n_features": n_f,
            "n_per_group": n_s,
            "n_total": 2 * n_s,
            "stage": "TOTAL",
            "t_s": round(total, 3),
        })

    return pd.DataFrame(records)


def main() -> None:
    print("PhenoRewire -- Component-Level Performance Profile")
    print("=" * 60)
    if not HAS_PSUTIL:
        print("(psutil not installed -- install with: pip install psutil)")
    print(f"N_PERMUTATIONS={N_PERMUTATIONS}, N_SELECTED={N_SELECTED} (post-selection feature count)")
    print()

    df = run_profiles()

    print("\n\nSummary table (seconds):")
    print("-" * 60)
    pivot = df.pivot(index="dataset", columns="stage", values="t_s")
    stage_order = ["permutation_test", "network_build", "rewiring", "TOTAL"]
    pivot = pivot.reindex(columns=[s for s in stage_order if s in pivot.columns])
    print(pivot.to_string())

    out_path = Path(__file__).parent / "performance_profile_results.csv"
    df.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")

    # Bottleneck annotation
    total_rows = df[df["stage"] == "TOTAL"].sort_values("n_features")
    perm_rows  = df[df["stage"] == "permutation_test"].sort_values("n_features")

    print("\nBottleneck analysis:")
    if len(perm_rows) >= 2:
        t_small = perm_rows.iloc[0]["t_s"]
        t_large = perm_rows.iloc[-1]["t_s"]
        n_small = perm_rows.iloc[0]["n_features"]
        n_large = perm_rows.iloc[-1]["n_features"]
        if t_small > 0:
            print(f"  Permutation test scales {n_large/n_small:.1f}x features -> "
                  f"{t_large/t_small:.1f}x time (expected ~{n_large/n_small:.1f}x linear)")

    if len(total_rows) >= 1:
        largest = total_rows.iloc[-1]
        t_total = largest["t_s"]
        print(f"\n  Largest config ({largest['dataset']}): {t_total:.1f}s total")
        if t_total < 60:
            print("  -> Interactive range (< 1 min).")
        elif t_total < 300:
            print("  -> Fast batch range (1-5 min).")
        elif t_total < 600:
            print("  -> Moderate batch range (5-10 min).")
        else:
            print(f"  -> Long run ({t_total/60:.0f} min). For 1000+ features, increase "
                  "N_PERMUTATIONS_PHENO to 1000 for better p-value resolution.")

    print("\nPractical guidance:")
    print("  - For <= 200 features: any N_PERMUTATIONS >= 200 is sufficient")
    print("  - For 500+ features:   use N_PERMUTATIONS >= 1000 (BH-FDR minimum p = 0.001)")
    print("  - Network construction bottleneck is O(n_selected^2) — keep n_selected <= 50")


if __name__ == "__main__":
    main()
