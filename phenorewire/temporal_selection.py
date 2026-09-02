# ── PhenoRewire · Step 2b — Temporal feature selection ───────────────────────
# Selects features associated with a continuous time variable via Spearman
# correlation, empirical permutation test, and Benjamini-Hochberg FDR correction.
# Reports directionality as higher abundance at one timepoint vs the other.
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import logging

from .utils import save_json
from .stats import adaptive_fdr_selection, spearman_permutation_test

logger = logging.getLogger(__name__)


def temporal_selection(
    X: pd.DataFrame,
    time: np.ndarray,
    feat_anno: pd.DataFrame,
    outdir: Path,
    n_permutations: int,
    fdr_alpha: float,
    adaptive_fdr: bool,
    min_selected_features: int,
    fdr_alpha_ceiling: float,
    seed: int = 42,
    min_features_hard_stop: int = 5,
) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info("Running temporal correlation selection...")

    feats = X.index.astype(str).tolist()
    mat = X.values.astype(float)
    n_perm = int(n_permutations)
    if n_perm < 10:
        raise ValueError("n_permutations must be >= 10 for stability.")

    r_obs, p_emp = spearman_permutation_test(
        mat,
        np.asarray(time, dtype=float),
        n_permutations=n_perm,
        seed=int(seed),
    )

    _, qvals, *_ = multipletests(p_emp, alpha=float(fdr_alpha_ceiling), method="fdr_bh")
    rej, selection_meta = adaptive_fdr_selection(
        qvals,
        base_alpha=float(fdr_alpha),
        adaptive=bool(adaptive_fdr),
        min_selected=int(min_selected_features),
        alpha_ceiling=float(fdr_alpha_ceiling),
    )

    n_selected = int(selection_meta["selected_feature_count"])
    if n_selected < int(min_features_hard_stop):
        adaptive_note = (
            f"  Adaptive FDR relaxation was {'enabled' if adaptive_fdr else 'disabled'}; "
            f"ceiling alpha = {fdr_alpha_ceiling}.\n"
            if adaptive_fdr else
            "  ADAPTIVE_SELECTION_THRESHOLDS is False; enable it or increase sample size.\n"
        )
        raise ValueError(
            f"Temporal feature selection hard stop: only {n_selected} feature(s) passed FDR "
            f"threshold (q <= {selection_meta['applied_fdr_alpha']:.4f}), fewer than "
            f"MIN_FEATURES_HARD_STOP = {min_features_hard_stop}.\n"
            + adaptive_note +
            "  Consider: more time points, fewer features, a higher FDR alpha, or enabling "
            "ADAPTIVE_SELECTION_THRESHOLDS."
        )

    res = pd.DataFrame({
        "feature_id": feats,
        "r_time": r_obs,
        "p_empirical": p_emp,
        "q_time": qvals,
        "significant": rej,
    }).set_index("feature_id")

    res_annot = res.join(feat_anno, how="left")
    res_annot.to_csv(outdir / "temporal_features_all.csv")

    selected = res_annot[res_annot["significant"]].copy()
    selected.to_csv(outdir / "temporal_features_selected.csv")

    summary = {
        "n_features_total": int(len(feats)),
        "n_features_selected": int(selected.shape[0]),
        "fdr_alpha": float(fdr_alpha),
        "applied_fdr_alpha": float(selection_meta["applied_fdr_alpha"]),
        "adaptive_fdr_used": bool(selection_meta["adaptive_fdr_used"]),
        "min_selected_target": int(selection_meta["min_selected_target"]),
        "selected_at_base_alpha": int(selection_meta["selected_at_base_alpha"]),
        "n_permutations": int(n_perm),
    }
    save_json(summary, outdir / "temporal_selection_summary.json")

    logger.info(f"Temporal selection done. Selected {selected.shape[0]} features.")
    return res_annot, selected, summary
