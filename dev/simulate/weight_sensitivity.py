"""
Triage weight sensitivity analysis for PhenoRewire.

Sweeps TRIAGE_TOPOLOGY_WEIGHT / TRIAGE_SELECTION_WEIGHT across [0.1/0.9 … 0.9/0.1]
using a synthetic dataset engineered so that topology signal and selection signal
rank features *differently*.

Design
------
* N_FEATURES = 25; N_SAMPLES_PER_GROUP = 30
* Ground truth = F000–F007 (8 features)
  - F000–F003 (GT_HUB):   dense hub block (|r| ~= 0.62) + moderate phenotype shift (+0.8)
                           -> high topology, moderate selection support
  - F004–F007 (GT_SEL):   no hub structure + strong phenotype shift (+1.4)
                           -> lower topology, high selection support
* F008–F012 (BG_HUB):     same hub block as F000–F003 + weak phenotype shift (+0.6)
                           -> high topology, lower selection support (non-GT, DECOY)
* F013–F024 (QUIET_BG):   no structure, no shift -> not selected, absent from triage

Expected AUROC behaviour
------------------------
* topology_weight=0.1 (selection-heavy): GT_SEL and GT_HUB both rank high due to strong
  association; BG_HUB (weak association) ranks below most GT -> AUROC ~= 0.8
* topology_weight=0.9 (topology-heavy): hub features (GT_HUB + BG_HUB) score high;
  GT_SEL (no hub) score low -> some BG_HUB features outrank GT_SEL -> AUROC ~= 0.6
* AUROC range >= 0.10 across weight pairs (crossover at topology_weight ~= 0.22)

Usage::

    python -m tests.simulate.weight_sensitivity        # from repo root
    python tests/simulate/weight_sensitivity.py        # direct

Output:
    - prints a summary table
    - saves weight_sensitivity_results.csv next to this script
    - saves weight_sensitivity_auroc.png if matplotlib is available
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from phenorewire.stats import fast_spearman_permutation_test
from phenorewire.networks import correlation_network, compute_network_metrics
from phenorewire.triage import export_network_triage

import tempfile
import os

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

N_HUB_GT = 4       # F000–F003: hub structure + moderate selection
N_SEL_GT = 4       # F004–F007: no hub + strong selection
N_HUB_BG = 5       # F008–F012: hub structure + weak selection (decoy)
N_QUIET_BG = 12    # F013–F024: no structure, no selection
N_FEATURES = N_HUB_GT + N_SEL_GT + N_HUB_BG + N_QUIET_BG  # 25
N_IMPORTANT = N_HUB_GT + N_SEL_GT                          # 8
N_SAMPLES_PER_GROUP = 30
SEED = 42

WEIGHT_PAIRS = [
    (0.1, 0.9),
    (0.3, 0.7),
    (0.5, 0.5),
    (0.7, 0.3),
    (0.9, 0.1),
]

# Hub block indices: F000-F003 + F008-F012
_HUB_IDX = list(range(N_HUB_GT)) + list(range(N_HUB_GT + N_SEL_GT, N_HUB_GT + N_SEL_GT + N_HUB_BG))


def _build_corr_matrix(n: int) -> np.ndarray:
    """
    Correlation matrix: a single 9-node hub block (indices _HUB_IDX) at r=0.62,
    all other pairs at r=0 (uncorrelated).

    Target |r| range: 0.55–0.70 within hub (signal), 0.10–0.25 off-hub (background).
    Off-hub is set to 0 here so the background never dominates topology.
    """
    C = np.eye(n)
    for i in _HUB_IDX:
        for j in _HUB_IDX:
            if i != j:
                C[i, j] = 0.62
    # Diagonal regularisation to guarantee PD
    eigvals = np.linalg.eigvalsh(C)
    if eigvals.min() <= 0:
        C += np.eye(n) * (-eigvals.min() + 0.01)
    return C


def simulate_dataset(seed: int = SEED) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """
    Generate a synthetic (features × samples) DataFrame with two groups.

    Returns
    -------
    X : pd.DataFrame, shape (N_FEATURES, 2 * N_SAMPLES_PER_GROUP)
    y : np.ndarray of int, binary group label
    ground_truth_ids : list[str]
    """
    rng = np.random.default_rng(seed)
    n_total = 2 * N_SAMPLES_PER_GROUP

    C = _build_corr_matrix(N_FEATURES)
    L = np.linalg.cholesky(C)

    # Draw correlated samples (both groups share the same covariance)
    Z = rng.standard_normal((n_total, N_FEATURES)) @ L.T  # (n_total, N_FEATURES)

    # Mean shifts in CASE group only (indices N_SAMPLES_PER_GROUP onward)
    shift = np.zeros(N_FEATURES)
    shift[:N_HUB_GT]                                         = 0.8   # F000-F003
    shift[N_HUB_GT : N_HUB_GT + N_SEL_GT]                   = 1.4   # F004-F007
    shift[N_HUB_GT + N_SEL_GT : N_HUB_GT + N_SEL_GT + N_HUB_BG] = 0.6   # F008-F012
    # F013-F024: 0 (no shift)
    Z[N_SAMPLES_PER_GROUP:] += shift

    X = Z.T  # (N_FEATURES, n_total)
    feature_ids = [f"F{i:03d}" for i in range(N_FEATURES)]
    sample_ids = [f"S{i:03d}" for i in range(n_total)]
    X_df = pd.DataFrame(X, index=feature_ids, columns=sample_ids)
    y = np.array([0] * N_SAMPLES_PER_GROUP + [1] * N_SAMPLES_PER_GROUP, dtype=int)
    ground_truth_ids = feature_ids[:N_IMPORTANT]
    return X_df, y, ground_truth_ids


def compute_triage_scores(
    X: pd.DataFrame,
    y: np.ndarray,
    *,
    topology_weight: float,
    selection_weight: float,
    seed: int = SEED,
    outdir: Path | None = None,
) -> pd.DataFrame:
    """
    Run feature selection + network + triage on a synthetic dataset.

    Returns a DataFrame with feature_id and triage_score columns.
    """
    from statsmodels.stats.multitest import multipletests
    from phenorewire.stats import adaptive_fdr_selection

    mat = X.values.astype(float)

    # Feature selection: permutation Spearman
    r_obs, p_emp = fast_spearman_permutation_test(mat, y, n_permutations=200, seed=seed)
    _, qvals, _, _ = multipletests(p_emp, alpha=0.2, method="fdr_bh")

    rej, sel_meta = adaptive_fdr_selection(
        qvals, base_alpha=0.2, adaptive=True, min_selected=5, alpha_ceiling=0.3
    )

    selected_ids = X.index[rej].tolist()
    if not selected_ids:
        order = np.argsort(-np.abs(np.where(np.isfinite(r_obs), r_obs, 0.0)))
        selected_ids = X.index[order[:10]].tolist()

    selected_feat = pd.DataFrame(
        {"q_pheno": qvals[rej]},
        index=pd.Index(selected_ids, name="feature_id"),
    )

    # Network construction
    X_sel = X.loc[selected_ids]
    _, edges_df, _ = correlation_network(
        X_sel,
        abs_r_thr=0.4,
        fdr_alpha=0.2,
        adaptive_thresholds=True,
        min_edges=5,
        max_density=0.5,
        floor_abs_r_thr=0.3,
        ceil_fdr_alpha=0.3,
    )
    nodes_df, community_df, network_summary_df = compute_network_metrics(selected_feat, edges_df, seed=seed)

    if outdir is None:
        tmp = tempfile.mkdtemp()
        outdir = Path(tmp)

    export_network_triage(
        nodes_df,
        community_df,
        network_summary_df,
        outdir,
        top_n=N_FEATURES,
        topology_weight=topology_weight,
        selection_weight=selection_weight,
    )

    triage_path = outdir / "triage" / "node_priority_ranking.csv"
    if triage_path.exists():
        df = pd.read_csv(triage_path)
        if "feature_id" not in df.columns and df.columns[0] != "feature_id":
            df = df.rename(columns={df.columns[0]: "feature_id"})
        return df[["feature_id", "triage_score"]].copy()

    return pd.DataFrame(columns=["feature_id", "triage_score"])


def run_sensitivity(seed: int = SEED) -> pd.DataFrame:
    """Sweep weight pairs and compute AUROC for ground-truth node recovery."""
    X, y, gt_ids = simulate_dataset(seed=seed)
    gt_set = set(gt_ids)

    records = []
    for t_w, s_w in WEIGHT_PAIRS:
        triage_df = compute_triage_scores(X, y, topology_weight=t_w, selection_weight=s_w, seed=seed)
        if triage_df.empty:
            records.append({"topology_weight": t_w, "selection_weight": s_w, "auroc": float("nan"), "n_features_ranked": 0})
            continue

        all_ids = triage_df["feature_id"].tolist()
        scores = triage_df["triage_score"].tolist()
        binary = [1 if fid in gt_set else 0 for fid in all_ids]
        if sum(binary) == 0 or sum(binary) == len(binary):
            auroc = float("nan")
        else:
            auroc = float(roc_auc_score(binary, scores))

        records.append({
            "topology_weight": t_w,
            "selection_weight": s_w,
            "auroc": round(auroc, 4),
            "n_features_ranked": len(all_ids),
        })

    return pd.DataFrame(records)


def main() -> None:
    print("PhenoRewire - Triage Weight Sensitivity Analysis")
    print("=" * 56)
    print(f"Synthetic dataset: {N_FEATURES} features, {N_SAMPLES_PER_GROUP} samples/group")
    print(f"Ground-truth important features: {N_IMPORTANT} (F000-F{N_IMPORTANT-1:03d})")
    print("Hub block: F000-F003 (GT) + F008-F012 (BG decoy) -> topology/selection tension built in")
    print()

    results = run_sensitivity(seed=SEED)
    print(results.to_string(index=False))

    out_path = Path(__file__).parent / "weight_sensitivity_results.csv"
    results.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")

    valid = results.dropna(subset=["auroc"])
    if valid.empty:
        print("\nAll AUROC values are NaN - check that GT and non-GT features both appear in triage output.")
        return

    best_row = valid.loc[valid["auroc"].idxmax()]
    worst_row = valid.loc[valid["auroc"].idxmin()]
    auroc_range = float(valid["auroc"].max() - valid["auroc"].min())

    print(f"\nBest:  topology={best_row['topology_weight']}, selection={best_row['selection_weight']}"
          f" -> AUROC = {best_row['auroc']:.4f}")
    print(f"Worst: topology={worst_row['topology_weight']}, selection={worst_row['selection_weight']}"
          f" -> AUROC = {worst_row['auroc']:.4f}")
    print(f"Range: {auroc_range:.4f}")

    if auroc_range >= 0.10:
        print(f"\nConclusion: Weights MATTER (range {auroc_range:.4f} >= 0.10). PASS")
        print(f"Best weight pair: topology={best_row['topology_weight']}, "
              f"selection={best_row['selection_weight']}.")
        print("Selection-heavy weights (0.1/0.9) favour strongly-associated features.")
        print("Topology-heavy weights (0.9/0.1) favour hub features regardless of q-value.")
        print("The 0.7/0.3 default balances both signals.")
    else:
        print(f"\nConclusion: AUROC range {auroc_range:.4f} < 0.10 - weight sensitivity FAIL.")
        print("Check simulation design: GT_SEL and BG_HUB topology/selection scores may overlap.")

    # Plot if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(
            [f"{r['topology_weight']}/{r['selection_weight']}" for _, r in results.iterrows()],
            results["auroc"],
            color=["#2196F3" if abs(r["topology_weight"] - 0.7) < 0.01 else "#90CAF9"
                   for _, r in results.iterrows()],
            edgecolor="white",
        )
        ax.axhline(0.5, color="red", linestyle="--", linewidth=0.8, label="random (0.5)")
        ax.set_xlabel("topology_weight / selection_weight")
        ax.set_ylabel("AUROC")
        ax.set_title("Triage weight sensitivity: AUROC for GT node recovery")
        ax.set_ylim(0, 1.05)
        ax.legend()
        plt.tight_layout()
        png_path = Path(__file__).parent / "weight_sensitivity_auroc.png"
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
        print(f"Plot saved to: {png_path}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
