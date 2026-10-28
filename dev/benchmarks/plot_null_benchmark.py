"""
Supplementary figure for null_benchmark.py.

Row 1: fraction of replicates in which each method claims the two states differ
in correlation structure.  In N and A the right answer is no (false-positive
rate); in D and B it is power.

  PhenoRewire 0.2.0       abundance + connectivity selection, rewiring null p <= 0.05
  abundance-only          selection by abundance only (as before 0.2.0) + rewiring null
  raw count               at least one state-specific pair, no null (the
                          decision PhenoRewire implied before 0.2.0)
  DiffCorr                at least one pair at BH q <= 0.05

Row 2: mean number of pairs reported as rewired that are not in the truth.
For the two PhenoRewire modes pairs count only in replicates whose comparison
passes the null, since otherwise nothing is reported.

Row 3 (D and B, where metabolites truly rewire): fraction of the rewired
metabolites each method selects.  Differential abundance (Mann-Whitney, BH) is
the reference: what a standard analysis would report.

Usage:  python dev/benchmarks/plot_null_benchmark.py [runs.csv] [--out stem]
Writes <stem>.svg, <stem>.png and <stem>.csv (the plotted values).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
TRUTH_PAIRS = 20 * 19 // 2
SCENARIOS = [
    ("N", "N  no difference", False),
    ("A", "A  abundance only", False),
    ("D", "D  abundance + rewiring", True),
    ("B", "B  rewiring only", True),
]
BLUE, AQUA, ORANGE = "#2a78d6", "#1baf7a", "#eb6834"   # reference slots 1-3, validated all-pairs
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
# key, label, colour, marker, linestyle, marker filled, x offset.  The two
# PhenoRewire modes share a hue (same tool) and differ in line and marker.
SERIES = [
    ("pr_union", "PhenoRewire 0.2.0", BLUE, "o", "-", True, 0.0),
    ("pr_abund", "abundance-only selection (before 0.2.0)", BLUE, "o", "--", False, -0.3),
    ("pr_raw", "raw state-specific count", AQUA, "s", ":", True, -0.6),
    ("diffcorr", "DiffCorr", ORANGE, "^", "-", True, 0.3),
    ("dgca", "DGCA", ORANGE, "v", "--", False, 0.45),
    ("da", "differential abundance", MUTED, "D", "-", True, 0.6),
]
ROWS = [("claims", "fraction of replicates\nclaiming rewiring"),
        ("false_pairs", "false pairs reported\n(mean)"),
        ("feature_recall", "rewired metabolites\nselected (fraction)")]


def tidy(runs: pd.DataFrame) -> pd.DataFrame:
    runs = runs[runs["status"] == "ok"].copy()
    hits = runs["pair_recall"].fillna(0.0) * np.where(runs["scenario"].isin(["D", "B"]), TRUTH_PAIRS, 0)
    runs["false_pairs"] = runs["pairs_called"].fillna(0) - hits

    parts = []
    for method, key in (("phenorewire_union", "pr_union"), ("phenorewire", "pr_abund")):
        pr = runs[runs["method"] == method]
        sig = pr["comparison_significant"].astype(bool)
        parts.append(pr.assign(series=key, claims=sig.astype(float),
                               false_pairs=pr["false_pairs"].where(sig, 0.0),
                               feature_recall=pr["selected_recall"]))
    pr = runs[runs["method"] == "phenorewire"]
    parts.append(pr.assign(series="pr_raw", claims=(pr["pairs_called"] > 0).astype(float),
                           feature_recall=np.nan))
    dc = runs[runs["method"] == "diffcorr"]
    parts.append(dc.assign(series="diffcorr", claims=(dc["pairs_called"] > 0).astype(float), feature_recall=np.nan))
    dg = runs[runs["method"] == "dgca"]
    parts.append(dg.assign(series="dgca", claims=(dg["pairs_called"] > 0).astype(float), feature_recall=np.nan))
    da = runs[runs["method"] == "da"]
    parts.append(da.assign(series="da", claims=np.nan, false_pairs=np.nan, feature_recall=da["selected_recall"]))
    long = pd.concat(parts)
    return (long.groupby(["series", "scenario", "n_per_state"])
            .agg(claims=("claims", "mean"), false_pairs=("false_pairs", "mean"),
                 feature_recall=("feature_recall", "mean"), reps=("rep", "count"))
            .reset_index())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="?", type=Path, default=HERE / "null_benchmark_runs.csv")
    ap.add_argument("--out", type=Path, default=HERE / "null_benchmark_figure")
    args = ap.parse_args()

    agg = tidy(pd.read_csv(args.runs))
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                         "svg.fonttype": "none"})
    fig, axes = plt.subplots(3, 4, figsize=(7.2, 5.6), sharex=True, sharey="row")
    ymax_false = max(1.0, np.nanmax(agg["false_pairs"].to_numpy()))
    n_lo, n_hi = agg["n_per_state"].min(), agg["n_per_state"].max()

    for col, (sc, title, rewired) in enumerate(SCENARIOS):
        for row, (metric, ylabel) in enumerate(ROWS):
            ax = axes[row, col]
            ax.grid(axis="y", color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.set_xlim(n_lo - 1.5, n_hi + 1.5)
            ax.axvspan(n_lo - 1.5, 10, color=GRID, alpha=0.45, linewidth=0, zorder=0)
            if row == 0:
                ax.set_title(title, fontsize=8, color=INK, loc="left", pad=12)
                ax.text(0.0, 1.015, "claims here = power" if rewired else "claims here = false positives",
                        transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5, color=MUTED, style="italic")
                if not rewired:
                    ax.axhline(0.05, color=MUTED, linewidth=0.8, linestyle=":")
            if row == 2 and not rewired:
                ax.text(0.5, 0.5, "no rewired\nmetabolites", transform=ax.transAxes, ha="center",
                        va="center", fontsize=7, color=MUTED)
            else:
                for key, label, color, marker, ls, filled, dx in SERIES:
                    d = agg[(agg["series"] == key) & (agg["scenario"] == sc)].sort_values("n_per_state")
                    d = d[np.isfinite(d[metric].astype(float))]
                    if d.empty:
                        continue
                    ax.plot(d["n_per_state"] + dx, d[metric], color=color, linestyle=ls, linewidth=1.8,
                            marker=marker, markersize=4.5, markeredgewidth=1.0,
                            markerfacecolor=color if filled else "white",
                            markeredgecolor=color if not filled else "white", zorder=3)
            if row == 1:
                ax.set_yscale("symlog", linthresh=1)
                ax.set_ylim(0, ymax_false * 1.5)
            else:
                ax.set_ylim(-0.03, 1.03)
            if row == 2:
                ax.set_xlabel("samples per state", color=INK)
            if col == 0:
                ax.set_ylabel(ylabel)

    axes[2, 0].set_xticks([n for n in sorted(agg["n_per_state"].unique()) if n != 6])
    handles = [plt.Line2D([], [], color=c, linestyle=ls, linewidth=1.8, marker=m, markersize=4.5,
                          markerfacecolor=c if f else "white", markeredgecolor=c if not f else "white")
               for _, _, c, m, ls, f, _ in SERIES]
    labels = [label for _, label, *_ in SERIES]
    handles.append(matplotlib.patches.Patch(color=GRID, alpha=0.8))
    labels.append("< 10 per state")
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=7,
               bbox_to_anchor=(0.5, -0.005), labelcolor=INK, handlelength=2.6, columnspacing=1.2)
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    for ext in ("svg", "png"):
        fig.savefig(args.out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    agg.to_csv(args.out.with_suffix(".csv"), index=False, float_format="%.3f")
    print(f"wrote {args.out}.svg/.png/.csv")


if __name__ == "__main__":
    main()
