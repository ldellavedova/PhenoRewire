# ── PhenoRewire · Step 2b — Temporal feature selection ───────────────────────
# Selects features associated with a continuous time variable via Spearman
# correlation, empirical permutation test, and Benjamini-Hochberg FDR correction.
# Directionality (which timepoint a feature is higher at) is added by run.py.
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import logging
import numpy as np
import pandas as pd

from .selection import run_permutation_selection

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
    """
    Time-driven feature selection via Spearman correlation with permutation-based empirical p-values.

    Assumptions:
      - X is features x samples (index = feature_id, columns = sample IDs)
      - time is a numeric 1D array aligned to X columns order
      - X is already preprocessed upstream

    Outputs (in outdir):
      - temporal_features_all.csv
      - temporal_features_selected.csv
      - temporal_selection_summary.json
    """
    logger.info("Running temporal correlation selection...")

    time_arr = np.asarray(time, dtype=float)
    if time_arr.ndim != 1:
        raise ValueError("time must be a 1D array.")
    if time_arr.shape[0] != X.shape[1]:
        raise ValueError(
            f"time length ({time_arr.shape[0]}) must match number of samples in X ({X.shape[1]})."
        )

    return run_permutation_selection(
        X,
        time_arr,
        feat_anno,
        outdir,
        mode="Temporal",
        r_col="r_time",
        q_col="q_time",
        all_filename="temporal_features_all.csv",
        selected_filename="temporal_features_selected.csv",
        summary_filename="temporal_selection_summary.json",
        n_permutations=n_permutations,
        fdr_alpha=fdr_alpha,
        adaptive_fdr=adaptive_fdr,
        min_selected_features=min_selected_features,
        fdr_alpha_ceiling=fdr_alpha_ceiling,
        seed=seed,
        min_features_hard_stop=min_features_hard_stop,
        remedy="more time points, fewer features",
    )
