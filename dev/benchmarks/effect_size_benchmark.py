"""How large must a reorganization be for PhenoRewire to see it at small sample size?

The main benchmark simulates one effect: a block of 20 of 200 metabolites (10%)
that is correlated at rho = 0.85 in one state and uncorrelated in the other. That
answers whether the framework is safe, not what it can still resolve when an
experiment has five or six samples per state, which is the regime of a
strain-dropout design.

Here the same scenario is run with the effect growing along two axes, how much of
the metabolome is involved and how tightly it is coupled, at the sample sizes
below the operating range. Marginal distributions are held identical between the
two states throughout (scenario B), so abundance testing has nothing to find and
only the reorganization can be detected.

    python dev/benchmarks/effect_size_benchmark.py [--reps 20]

Output: effect_size_runs.csv, effect_size_summary.csv
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import null_benchmark as nb  # noqa: E402

P_FEATURES = nb.P_FEATURES


def simulate(n: int, block: int, rho: float, shift: float, rng):
    """A block that loses its correlation and, if shift > 0, also changes in level.

    shift = 0 is pure rewiring with identical marginals (scenario B of the main
    benchmark); shift > 0 is the case the strain-dropout comparisons of the study
    fall into, where the reorganized metabolites also differ in abundance and can
    therefore be selected below 14 samples per state.
    """
    ref = rng.standard_normal((n, P_FEATURES))
    case = rng.standard_normal((n, P_FEATURES))
    idx = np.arange(block)
    if block > 0:
        ref[:, idx] = nb._block(n, block, rho, rng)
        case[:, idx] = nb._block(n, block, rho if shift < 0 else 0.0, rng)
    if shift == 0:
        for arr in (ref, case):                  # identical marginals
            sub = arr[:, idx]
            arr[:, idx] = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=1)
    else:
        case[:, idx] += abs(shift)
    truth = {(int(i), int(j)) for i in idx for j in idx if i < j} if (block > 0 and shift >= 0) else set()
    return np.exp(ref.T), np.exp(case.T), truth, set(idx.tolist())


def run_pr(a, b, *, sel_alpha, edge_alpha, abs_r, perm, null_perm, seed):
    """PhenoRewire with selection and network thresholds set independently.

    The published benchmark runs every method at one matched false-discovery
    rate, which is the right way to compare methods. The runs of this study use
    two levels, as the configuration files record: selection at its own ceiling,
    edges at |rho| and their own FDR. This function reproduces that.
    """
    import numpy as np
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    from phenorewire.networks import compute_rewiring, correlation_network
    from phenorewire.connectivity import CONNECTIVITY_MIN_N, connectivity_selection
    from phenorewire.nulls import rewiring_label_null
    from phenorewire.stats import spearman_permutation_test

    X = np.log2(np.hstack([a, b]) + 1.0)
    y = np.r_[np.zeros(a.shape[1]), np.ones(b.shape[1])]
    _, p = spearman_permutation_test(X, y, n_permutations=perm, seed=seed)
    q = multipletests(p, method="fdr_bh")[1]
    picked_set = set(np.flatnonzero(q <= sel_alpha).tolist())
    if min(a.shape[1], b.shape[1]) >= CONNECTIVITY_MIN_N:
        fe = [f"F{i:03d}" for i in range(X.shape[0])]
        tab, _ = connectivity_selection(
            pd.DataFrame(X[:, : a.shape[1]], index=fe), pd.DataFrame(X[:, a.shape[1]:], index=fe),
            top_k=10, n_permutations=200, fdr_alpha=0.05, seed=seed)
        picked_set |= set(np.flatnonzero(tab["selected_connectivity"].to_numpy()).tolist())
    picked = np.array(sorted(picked_set), dtype=int)
    out = {"n_selected": int(picked.size), "pairs": set(), "null_p": float("nan"),
           "obs_over_null": float("nan")}
    if picked.size < 3:
        return out
    feats = [f"F{i:03d}" for i in picked]
    Xa = pd.DataFrame(X[picked, : a.shape[1]], index=feats)
    Xb = pd.DataFrame(X[picked, a.shape[1]:], index=feats)
    kw = dict(abs_r_thr=abs_r, fdr_alpha=edge_alpha, adaptive_thresholds=False)
    _, ea, _ = correlation_network(Xa, **kw)
    _, eb, _ = correlation_network(Xb, **kw)
    _, only_a, only_b, _ = compute_rewiring(ea, eb, label_a="a", label_b="b")
    for df in (only_a, only_b):
        for u, v in zip(df["node_a"], df["node_b"]):
            i, j = int(u[1:]), int(v[1:])
            out["pairs"].add((min(i, j), max(i, j)))
    summary, _ = rewiring_label_null(Xa, Xb, thresholds_a=(abs_r, edge_alpha),
                                     thresholds_b=(abs_r, edge_alpha),
                                     n_permutations=null_perm, seed=seed)
    out["null_p"] = summary["p_empirical"]
    out["obs_over_null"] = summary["obs_over_null"]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--n-grid", type=int, nargs="+", default=[5, 6, 10, 14])
    ap.add_argument("--blocks", type=int, nargs="+", default=[20, 40, 80])
    ap.add_argument("--rhos", type=float, nargs="+", default=[0.85, 0.95])
    ap.add_argument("--shifts", type=float, nargs="+", default=[0.0, 1.2, 2.0])
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--null-perm", type=int, default=200)
    ap.add_argument("--sel-alpha", type=float, default=0.20,
                    help="feature-selection FDR; 0.20 is the ceiling the strain-dropout runs used")
    ap.add_argument("--edge-alpha", type=float, default=0.05, help="edge FDR, as in the runs")
    ap.add_argument("--abs-r", type=float, default=0.70, help="|rho| for edges, as in the runs")
    ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--rscript", default=None, help='e.g. "conda run -n rbench Rscript"; adds DiffCorr and DGCA')
    args = ap.parse_args()

    rows = []
    t0 = time.perf_counter()
    for n in args.n_grid:
        for block in args.blocks:
            for rho in args.rhos:
              for shift in args.shifts:
                for rep in range(args.reps):
                    seed = args.seed + 7919 * rep + 31 * n + 13 * block + int(rho * 100) + int(shift * 10)
                    rng = np.random.default_rng(seed)
                    a, b, truth, truth_feats = simulate(n, block, rho, shift, rng)
                    pr = run_pr(a, b, sel_alpha=args.sel_alpha, edge_alpha=args.edge_alpha,
                                abs_r=args.abs_r, perm=args.perm, null_perm=args.null_perm, seed=seed)
                    comp = {}
                    if args.rscript:
                        import shlex
                        import os
                        rs = [t.strip('"') for t in shlex.split(args.rscript, posix=(os.name != "nt"))]
                        for name, fn in (("diffcorr", lambda: nb.run_diffcorr(a, b, 0.05, rs)),
                                         ("dgca", lambda: nb.run_dgca(a, b, 0.05, rs, 0))):
                            pairs, _ = fn()
                            comp[name] = bool(pairs) if pairs is not None else None
                    rec, prec = nb._recall_precision(pr["pairs"], truth)
                    rows.append({
                        "n_per_state": n, "block": block, "block_fraction": block / P_FEATURES,
                        "rho": rho, "shift": shift, "rep": rep,
                        "detected": bool(pr["null_p"] <= 0.05) if pr["null_p"] == pr["null_p"] else False,
                        "obs_over_null": pr["obs_over_null"], "null_p": pr["null_p"],
                        "n_selected": pr["n_selected"], "pair_recall": rec, "pair_precision": prec,
                        "diffcorr_claims": comp.get("diffcorr"), "dgca_claims": comp.get("dgca"),
                    })
                print(f"n={n:3d} block={block:3d} rho={rho:.2f} shift={shift:.1f}  "
                      f"detected {np.mean([r['detected'] for r in rows[-args.reps:]]):.0%}  "
                      f"({time.perf_counter()-t0:.0f}s)", flush=True)

    runs = pd.DataFrame(rows)
    runs.to_csv(HERE / "effect_size_runs.csv", index=False)
    summary = (runs.groupby(["n_per_state", "block", "block_fraction", "rho", "shift"])
               .agg(reps=("rep", "count"),
                    detected=("detected", "mean"),
                    diffcorr_claims=("diffcorr_claims", "mean"),
                    dgca_claims=("dgca_claims", "mean"),
                    median_obs_over_null=("obs_over_null", "median"),
                    median_selected=("n_selected", "median"),
                    pair_recall=("pair_recall", "mean"))
               .reset_index())
    summary.to_csv(HERE / "effect_size_summary.csv", index=False)
    print("\n", summary.to_string(index=False))


if __name__ == "__main__":
    main()
