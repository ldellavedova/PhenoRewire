"""
Ground-truth benchmark for PhenoRewire 0.2.0, including its rewiring null.

Every PhenoRewire number below is produced by the package's own functions
(stats.spearman_permutation_test, networks.correlation_network,
networks.compute_rewiring, nulls.rewiring_label_null), not by a
re-implementation, so the benchmark measures the released code.

Simulated two-state data, 200 features, one block of 20 carrying the signal:

  N  null          nothing differs between the states
  A  abundance     the block shifts its mean; its correlation is unchanged
  D  both          the block shifts its mean AND loses its correlation
                   (the ordinary biological case; what PhenoRewire targets)
  B  rewiring      the block loses its correlation; marginals identical

Three levels are scored:

  comparison  does the rewiring null call the comparison significant
              (p <= 0.05)?  In N and A that is a false positive; in D and B,
              power.  This is the decision a user acts on.
  pair        recall / precision of the raw state-specific edges against the
              truth (pairs inside the rewired block), and how many are called
              when the truth is zero.
  feature     recall / precision of features whose state-specific degree
              exceeds the null at BH q <= 0.05.

DiffCorr (Fukushima 2013), a published differential-correlation method, is
run through Rscript when available (--rscript), on identical data, scored at
the pair level at matched FDR.

Usage:
    python dev/benchmarks/null_benchmark.py                       # PhenoRewire only
    python dev/benchmarks/null_benchmark.py --rscript "conda run -n rbench Rscript"

Outputs (next to this script, or --outdir):
    null_benchmark_runs.csv     one row per method x scenario x n x replicate
    null_benchmark_summary.csv  means per method x scenario x n
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from phenorewire import __version__
from phenorewire.networks import compute_rewiring, correlation_network
from phenorewire.connectivity import CONNECTIVITY_MIN_N, connectivity_selection
from phenorewire.nulls import rewiring_label_null
from phenorewire.stats import spearman_permutation_test

P_FEATURES = 200
BLOCK = 20
MEAN_SHIFT = 1.2
RHO = 0.85
SCENARIOS = {
    "N": "no difference",
    "A": "mean shift only",
    "D": "mean shift and loss of correlation, same block",
    "B": "loss of correlation only, marginals identical",
}
COMPARATORS = Path(__file__).resolve().parent / "comparators"
ADAPTER = COMPARATORS / "diffcorr_adapter.R"          # kept: used by run_diffcorr
ADAPTER_DGCA = COMPARATORS / "dgca_adapter.R"


# ---------------------------------------------------------------- simulation
def _block(n: int, k: int, rho: float, rng) -> np.ndarray:
    """n samples x k features, equicorrelated at rho."""
    if rho <= 0:
        return rng.standard_normal((n, k))
    return np.sqrt(rho) * rng.standard_normal((n, 1)) + np.sqrt(1 - rho) * rng.standard_normal((n, k))


def simulate(scenario: str, n: int, rng):
    """Two states as features x samples (positive, log-normal-like), plus the truth."""
    ref = rng.standard_normal((n, P_FEATURES))
    case = rng.standard_normal((n, P_FEATURES))
    idx = np.arange(BLOCK)
    rewired = scenario in ("D", "B")

    if scenario in ("A", "D", "B"):
        ref[:, idx] = _block(n, BLOCK, RHO, rng)
        case[:, idx] = _block(n, BLOCK, 0.0 if rewired else RHO, rng)
    if scenario in ("A", "D"):
        case[:, idx] += MEAN_SHIFT
    if scenario == "B":
        for arr in (ref, case):
            sub = arr[:, idx]
            arr[:, idx] = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=1)

    truth = {(int(i), int(j)) for i in idx for j in idx if i < j} if rewired else set()
    return np.exp(ref.T), np.exp(case.T), truth, set(idx.tolist()) if rewired else set()


# --------------------------------------------------------------- PhenoRewire
def run_phenorewire(a: np.ndarray, b: np.ndarray, args, seed: int, mode: str = "abundance") -> dict:
    """Selection -> per-state networks -> rewiring -> label-permutation null.

    ``mode`` "abundance" selects by abundance only (PhenoRewire before 0.2.0);
    "union" adds connectivity-based selection from
    CONNECTIVITY_MIN_N samples per state, as v0.2.0 does automatically.
    """
    X = np.log2(np.hstack([a, b]) + 1.0)
    y = np.r_[np.zeros(a.shape[1]), np.ones(b.shape[1])]
    _, p = spearman_permutation_test(X, y, n_permutations=args.perm, seed=seed)
    q = multipletests(p, method="fdr_bh")[1]
    abundance = set(np.flatnonzero(q <= args.alpha).tolist())
    connectivity: set = set()
    if mode == "union" and min(a.shape[1], b.shape[1]) >= CONNECTIVITY_MIN_N:
        fe = [f"F{i:03d}" for i in range(X.shape[0])]
        tab, _ = connectivity_selection(
            pd.DataFrame(X[:, : a.shape[1]], index=fe), pd.DataFrame(X[:, a.shape[1]:], index=fe),
            top_k=10, n_permutations=200, fdr_alpha=args.alpha, seed=seed,
        )
        connectivity = set(np.flatnonzero(tab["selected_connectivity"].to_numpy()).tolist())
    picked = np.array(sorted(abundance | connectivity), dtype=int)
    out = {"n_selected": int(picked.size), "selected": set(picked.tolist()),
           "connectivity_only": connectivity - abundance, "pairs": set(), "null_p": np.nan,
           "obs_over_null": np.nan, "features_q05": set()}
    if picked.size < 3:
        return out

    feats = [f"F{i:03d}" for i in picked]
    Xa = pd.DataFrame(X[picked, : a.shape[1]], index=feats)
    Xb = pd.DataFrame(X[picked, a.shape[1]:], index=feats)
    kw = dict(abs_r_thr=args.abs_r, fdr_alpha=args.alpha, adaptive_thresholds=False)
    _, ea, _ = correlation_network(Xa, **kw)
    _, eb, _ = correlation_network(Xb, **kw)
    _, only_a, only_b, _ = compute_rewiring(ea, eb, label_a="a", label_b="b")
    for df in (only_a, only_b):
        for u, v in zip(df["node_a"], df["node_b"]):
            i, j = int(u[1:]), int(v[1:])
            out["pairs"].add((min(i, j), max(i, j)))

    summary, nodes = rewiring_label_null(
        Xa, Xb, thresholds_a=(args.abs_r, args.alpha), thresholds_b=(args.abs_r, args.alpha),
        n_permutations=args.null_perm, seed=seed,
    )
    out["null_p"] = summary["p_empirical"]
    out["obs_over_null"] = summary["obs_over_null"]
    out["features_q05"] = {int(f[1:]) for f in nodes.loc[nodes["q_bh"] <= 0.05, "feature_id"]}
    return out


# --------------------------------------------------- differential abundance
def run_da(a: np.ndarray, b: np.ndarray, alpha: float) -> set:
    """The standard per-feature test: Mann-Whitney U on log intensities, BH."""
    from scipy.stats import mannwhitneyu

    la, lb = np.log2(a + 1.0), np.log2(b + 1.0)
    p = np.array([mannwhitneyu(x, z, alternative="two-sided").pvalue for x, z in zip(la, lb)])
    return set(np.flatnonzero(multipletests(p, method="fdr_bh")[1] <= alpha).tolist())


# ----------------------------------------------------------------- DiffCorr
def run_r_pairs(adapter: Path, a: np.ndarray, b: np.ndarray, alpha: float, rscript,
                extra: list[str] | None = None) -> tuple[set | None, str]:
    """Run an R comparator that returns the feature pairs it calls."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "in.json").write_text(json.dumps({"a": a.tolist(), "b": b.tolist()}), encoding="utf-8")
        try:
            proc = subprocess.run(
                list(rscript) + [str(adapter), str(td / "in.json"), str(td / "out.json"), str(alpha)]
                + list(extra or []),
                capture_output=True, text=True, timeout=900,
            )
        except subprocess.TimeoutExpired:
            return None, "timeout"
        if proc.returncode != 0 or not (td / "out.json").is_file():
            tail = (proc.stderr or "").strip().splitlines()
            return None, tail[-1][:160] if tail else "failed"
        rows = json.loads((td / "out.json").read_text(encoding="utf-8")) or []
        return {(min(int(i), int(j)), max(int(i), int(j))) for i, j in rows}, ""


def run_diffcorr(a, b, alpha, rscript):
    return run_r_pairs(ADAPTER, a, b, alpha, rscript)


def run_dgca(a, b, alpha, rscript, nperm: int):
    return run_r_pairs(ADAPTER_DGCA, a, b, alpha, rscript, extra=[str(nperm)])


# ------------------------------------------------------------------ scoring
def _recall_precision(found: set, truth: set) -> tuple[float, float]:
    hit = len(found & truth)
    recall = hit / len(truth) if truth else np.nan
    precision = hit / len(found) if found else np.nan
    return recall, precision


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-grid", type=int, nargs="+", default=[5, 6, 10, 14, 20, 25])
    ap.add_argument("--scenarios", nargs="+", default=list(SCENARIOS))
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--perm", type=int, default=1000, help="selection permutations (published configs use 1000)")
    ap.add_argument("--null-perm", type=int, default=200, help="label permutations for the rewiring null")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--abs-r", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=20261014)
    ap.add_argument("--dgca-perm", type=int, default=0,
                    help="DGCA label permutations for its empirical FDR; 0 uses the BH-adjusted "
                         "asymptotic p, which matches how DiffCorr is scored")
    ap.add_argument("--rscript", default=None, help='e.g. "conda run -n rbench Rscript"; omit to skip DiffCorr')
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument("--resume", action="store_true", help="keep completed replicates from an earlier run")
    args = ap.parse_args()

    rscript = None
    if args.rscript:
        rscript = [t.strip('"') for t in shlex.split(args.rscript, posix=(os.name != "nt"))]
        probe = subprocess.run(rscript + ["-e", "suppressMessages(library(DiffCorr)); cat('ok')"],
                               capture_output=True, text=True)
        if probe.returncode != 0 or "ok" not in probe.stdout:
            print("DiffCorr not reachable; continuing with PhenoRewire only")
            rscript = None
    has_dgca = False
    if rscript:
        probe = subprocess.run(rscript + ["-e", "suppressMessages(library(DGCA)); cat('ok')"],
                               capture_output=True, text=True)
        has_dgca = probe.returncode == 0 and "ok" in probe.stdout

    print(f"PhenoRewire {__version__}; DiffCorr: {'yes' if rscript else 'no'}; "
          f"DGCA: {'yes' if has_dgca else 'no'} (nPerm={args.dgca_perm})")
    rows = []
    args.outdir.mkdir(parents=True, exist_ok=True)
    runs_path = args.outdir / "null_benchmark_runs.csv"
    methods = ["phenorewire", "phenorewire_union", "da"]
    if rscript:
        methods.append("diffcorr")
    if has_dgca:
        methods.append("dgca")
    done: dict = {}
    if args.resume and runs_path.is_file():
        # Every replicate is seeded from (scenario, n, rep) alone, so a resumed
        # run produces exactly the rows an uninterrupted one would have; only
        # the methods still missing for a replicate are computed.
        prev = pd.read_csv(runs_path)
        prev = prev[prev["status"] == "ok"]
        rows = prev.to_dict(orient="records")
        for key, g in prev.groupby(["scenario", "n_per_state", "rep"]):
            done[(str(key[0]), int(key[1]), int(key[2]))] = set(g["method"])
        print(f"resuming: {len(rows)} method-replicates already complete")
    for scenario in args.scenarios:
        for n in args.n_grid:
            for rep in range(args.reps):
                todo = [m for m in methods if m not in done.get((scenario, n, rep), set())]
                if not todo:
                    continue
                seed = args.seed + 7919 * rep + 31 * n + ord(scenario)
                rng = np.random.default_rng(seed)
                a, b, truth, truth_feats = simulate(scenario, n, rng)
                base = {"scenario": scenario, "n_per_state": n, "rep": rep}

                for method, mode in (("phenorewire", "abundance"), ("phenorewire_union", "union")):
                    if method not in todo:
                        continue
                    t0 = time.perf_counter()
                    pr = run_phenorewire(a, b, args, seed, mode=mode)
                    rec, prec = _recall_precision(pr["pairs"], truth)
                    frec, fprec = _recall_precision(pr["features_q05"], truth_feats)
                    srec, sprec = _recall_precision(pr["selected"], truth_feats)
                    significant = bool(pr["null_p"] <= 0.05) if pr["null_p"] == pr["null_p"] else False
                    rows.append({**base, "method": method,
                                 "n_selected": pr["n_selected"],
                                 "selected_recall": srec, "selected_precision": sprec,
                                 "connectivity_only_selected": len(pr["connectivity_only"]),
                                 "connectivity_only_true": len(pr["connectivity_only"] & truth_feats),
                                 "pairs_called": len(pr["pairs"]), "pair_recall": rec, "pair_precision": prec,
                                 "null_p": pr["null_p"], "obs_over_null": pr["obs_over_null"],
                                 "comparison_significant": significant,
                                 "features_called_q05": len(pr["features_q05"]),
                                 "feature_recall": frec, "feature_precision": fprec,
                                 "seconds": round(time.perf_counter() - t0, 3), "status": "ok"})

                if "da" in todo:
                    t0 = time.perf_counter()
                    sel = run_da(a, b, args.alpha)
                    srec, sprec = _recall_precision(sel, truth_feats)
                    rows.append({**base, "method": "da", "n_selected": len(sel),
                                 "selected_recall": srec, "selected_precision": sprec,
                                 "seconds": round(time.perf_counter() - t0, 3), "status": "ok"})

                if "dgca" in todo:
                    t0 = time.perf_counter()
                    pairs, note = run_dgca(a, b, args.alpha, rscript, args.dgca_perm)
                    if pairs is None:
                        rows.append({**base, "method": "dgca", "status": note})
                    else:
                        rec, prec = _recall_precision(pairs, truth)
                        rows.append({**base, "method": "dgca", "pairs_called": len(pairs),
                                     "pair_recall": rec, "pair_precision": prec,
                                     "comparison_significant": bool(pairs),
                                     "seconds": round(time.perf_counter() - t0, 3), "status": "ok"})

                if "diffcorr" in todo:
                    t0 = time.perf_counter()
                    pairs, note = run_diffcorr(a, b, args.alpha, rscript)
                    if pairs is None:
                        rows.append({**base, "method": "diffcorr", "status": note})
                    else:
                        rec, prec = _recall_precision(pairs, truth)
                        rows.append({**base, "method": "diffcorr", "pairs_called": len(pairs),
                                     "pair_recall": rec, "pair_precision": prec,
                                     "comparison_significant": bool(pairs),
                                     "seconds": round(time.perf_counter() - t0, 3), "status": "ok"})

            runs = pd.DataFrame(rows)
            runs.to_csv(runs_path, index=False)
            cur = runs[(runs["scenario"] == scenario) & (runs["n_per_state"] == n)]
            msg = ", ".join(
                f"{m}: sig={g['comparison_significant'].mean():.2f} pairs={g['pairs_called'].mean():.1f}"
                for m, g in cur.groupby("method") if m != "da"
            )
            print(f"  {scenario} n={n:3d}  {msg}", flush=True)

    runs = pd.DataFrame(rows)
    ok = runs[runs["status"] == "ok"]
    summary = (
        ok.groupby(["scenario", "n_per_state", "method"])
        .agg(reps=("rep", "count"),
             frac_comparisons_significant=("comparison_significant", "mean"),
             mean_pairs_called=("pairs_called", "mean"),
             pair_recall=("pair_recall", "mean"),
             pair_precision=("pair_precision", "mean"),
             mean_obs_over_null=("obs_over_null", "mean") if "obs_over_null" in ok else ("rep", "count"),
             feature_recall=("feature_recall", "mean") if "feature_recall" in ok else ("rep", "count"),
             feature_precision=("feature_precision", "mean") if "feature_precision" in ok else ("rep", "count"),
             n_selected=("n_selected", "mean"),
             selected_recall=("selected_recall", "mean"),
             selected_precision=("selected_precision", "mean"),
             connectivity_only_selected=("connectivity_only_selected", "mean"),
             connectivity_only_true=("connectivity_only_true", "mean"),
             median_seconds=("seconds", "median"))
        .reset_index()
    )
    summary.insert(1, "truth", summary["scenario"].map(SCENARIOS))
    summary.to_csv(args.outdir / "null_benchmark_summary.csv", index=False, float_format="%.3f")
    print(f"\nwrote {args.outdir / 'null_benchmark_summary.csv'}")


if __name__ == "__main__":
    main()
