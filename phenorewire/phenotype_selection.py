# ── PhenoRewire · Step 2a — Phenotype feature selection ──────────────────────
# Selects features associated with group identity via Spearman correlation,
# empirical permutation test, and Benjamini-Hochberg FDR correction.
# Reports group means and log2FC (case_vs_ref) for selected features.
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import logging
import numpy as np
import pandas as pd

from .selection import run_permutation_selection

logger = logging.getLogger(__name__)


def phenotype_selection(
    X: pd.DataFrame,
    y: np.ndarray,
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
    Phenotype-driven feature selection via Spearman correlation with permutation-based empirical p-values.
    Adds effect-size columns (mean_ref/mean_case/log2FC_case_vs_ref) for downstream pathway enrichment.

    Assumptions:
      - X is features x samples (index = feature_id, columns = sample IDs)
      - y is binary phenotype with REF=0, CASE=1 (aligned to X columns order)
      - X is already preprocessed (presence filter / normalization / optional log transform) upstream

    Outputs (in outdir):
      - phenotype_features_all.csv
      - phenotype_features_selected.csv
      - phenotype_selection_summary.json
    """
    logger.info("Running phenotype-driven selection...")

    mat = X.values.astype(float)

    y_arr = np.asarray(y).astype(int)
    if y_arr.ndim != 1:
        raise ValueError("y must be a 1D array.")
    if y_arr.shape[0] != mat.shape[1]:
        raise ValueError(
            f"y length ({y_arr.shape[0]}) must match number of samples in X ({mat.shape[1]})."
        )
    if not set(np.unique(y_arr)).issubset({0, 1}):
        raise ValueError("y must be binary encoded as 0/1 (REF=0, CASE=1).")

    mask_ref = y_arr == 0
    mask_case = y_arr == 1
    n_ref = int(mask_ref.sum())
    n_case = int(mask_case.sum())
    if n_ref < 2 or n_case < 2:
        raise ValueError(f"Need >=2 samples per group. Got n_ref={n_ref}, n_case={n_case}.")

    # If upstream did log2(x+1), then mean_case - mean_ref is already log2FC.
    mean_ref = np.nanmean(mat[:, mask_ref], axis=1)
    mean_case = np.nanmean(mat[:, mask_case], axis=1)

    return run_permutation_selection(
        X,
        y_arr,
        feat_anno,
        outdir,
        mode="Phenotype",
        r_col="r_pheno",
        q_col="q_pheno",
        all_filename="phenotype_features_all.csv",
        selected_filename="phenotype_features_selected.csv",
        summary_filename="phenotype_selection_summary.json",
        n_permutations=n_permutations,
        fdr_alpha=fdr_alpha,
        adaptive_fdr=adaptive_fdr,
        min_selected_features=min_selected_features,
        fdr_alpha_ceiling=fdr_alpha_ceiling,
        seed=seed,
        min_features_hard_stop=min_features_hard_stop,
        remedy="more samples, fewer features",
        extra_columns={
            "mean_ref": mean_ref.astype(float),
            "mean_case": mean_case.astype(float),
            "log2FC_case_vs_ref": (mean_case - mean_ref).astype(float),
        },
        extra_summary={"n_ref": n_ref, "n_case": n_case},
    )
