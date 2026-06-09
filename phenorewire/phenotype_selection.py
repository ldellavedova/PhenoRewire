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
from statsmodels.stats.multitest import multipletests

from .utils import save_json
from .stats import adaptive_fdr_selection, fast_spearman_permutation_test

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
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info("Running phenotype-driven selection...")

    # ----------------------------
    # Inputs / alignment checks
    # ----------------------------
    feats = X.index.astype(str).str.strip().tolist()
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

    # ----------------------------
    # Effect sizes
    # ----------------------------
    # If upstream did log2(x+1), then mean_case - mean_ref is already log2FC.
    mean_ref = np.nanmean(mat[:, mask_ref], axis=1)
    mean_case = np.nanmean(mat[:, mask_case], axis=1)
    log2fc_case_vs_ref = mean_case - mean_ref

    # ----------------------------
    # Permutation Spearman
    # ----------------------------
    n_perm = int(n_permutations)
    if n_perm < 10:
        raise ValueError("n_permutations must be >= 10 for stability.")
    r_obs, p_emp = fast_spearman_permutation_test(
        mat,
        y_arr,
        n_permutations=n_perm,
        seed=int(seed),
    )

    # ----------------------------
    # Multiple testing correction
    # ----------------------------
    # alpha here affects only the reject boolean array (unused).
    # q-values are re-filtered manually below against fdr_alpha.
    # Using fdr_alpha (not ceiling) makes the intent clear to reviewers.
    _, qvals, _, _ = multipletests(p_emp, alpha=float(fdr_alpha), method="fdr_bh")
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
            f"Phenotype feature selection hard stop: only {n_selected} feature(s) passed FDR "
            f"threshold (q <= {selection_meta['applied_fdr_alpha']:.4f}), fewer than "
            f"MIN_FEATURES_HARD_STOP = {min_features_hard_stop}.\n"
            + adaptive_note +
            "  Consider: more samples, fewer features, a higher FDR_ALPHA, or enabling "
            "ADAPTIVE_SELECTION_THRESHOLDS."
        )

    res = (
        pd.DataFrame(
            {
                "feature_id": feats,
                "r_pheno": r_obs.astype(float),
                "p_empirical": p_emp.astype(float),
                "q_pheno": qvals.astype(float),
                "significant": rej.astype(bool),
                "mean_ref": mean_ref.astype(float),
                "mean_case": mean_case.astype(float),
                "log2FC_case_vs_ref": log2fc_case_vs_ref.astype(float),
            }
        )
        .set_index("feature_id")
    )

    # ----------------------------
    # Robust annotation attachment + proper "name"
    # ----------------------------
    res.index = res.index.astype(str).str.strip()

    # Allow feat_anno to be None/empty without crashing
    if feat_anno is None:
        feat_anno2 = pd.DataFrame(index=res.index)
        feat_anno2.index.name = "feature_id"
    else:
        feat_anno2 = feat_anno.copy()

        # Ensure feat_anno2 is indexed by feature_id
        if feat_anno2.index.name != "feature_id":
            if "feature_id" in feat_anno2.columns:
                feat_anno2["feature_id"] = feat_anno2["feature_id"].astype(str).str.strip()
                feat_anno2 = feat_anno2.set_index("feature_id", drop=False)
            else:
                feat_anno2.index = feat_anno2.index.astype(str).str.strip()
                feat_anno2.index.name = "feature_id"

        # Normalize the name column from common aliases (priority order)
        if "name" not in feat_anno2.columns:
            for alt in [
                "NAME",
                "compound_name",
                "metabolite_name",
                "consensus_annotation",
                "annotation",
            ]:
                if alt in feat_anno2.columns:
                    feat_anno2 = feat_anno2.rename(columns={alt: "name"})
                    break

    res_annot = res.join(feat_anno2, how="left")

    # Guarantee a usable "name" column (prefer true name, fallback to feature_id)
    fallback = pd.Series(res_annot.index.astype(str), index=res_annot.index, name="name")
    if "name" not in res_annot.columns:
        res_annot["name"] = fallback
    else:
        # Clean + fallback, avoiding the literal string "nan"
        name_clean = res_annot["name"].replace({None: np.nan})
        name_clean = name_clean.where(name_clean.notna(), fallback)
        name_clean = name_clean.astype(str).replace({"nan": np.nan, "": np.nan})
        res_annot["name"] = name_clean.fillna(fallback)

    # ----------------------------
    # Save outputs
    # ----------------------------
    res_annot.to_csv(outdir / "phenotype_features_all.csv", index=True)

    selected = res_annot[res_annot["significant"]].copy()
    selected.to_csv(outdir / "phenotype_features_selected.csv", index=True)

    summary = {
        "n_features_total": int(len(feats)),
        "n_features_selected": int(selected.shape[0]),
        "fdr_alpha": float(fdr_alpha),
        "applied_fdr_alpha": float(selection_meta["applied_fdr_alpha"]),
        "adaptive_fdr_used": bool(selection_meta["adaptive_fdr_used"]),
        "min_selected_target": int(selection_meta["min_selected_target"]),
        "selected_at_base_alpha": int(selection_meta["selected_at_base_alpha"]),
        "n_permutations": int(n_perm),
        "n_ref": int(n_ref),
        "n_case": int(n_case),
    }
    save_json(summary, outdir / "phenotype_selection_summary.json")

    logger.info("Phenotype selection done. Selected %d features.", selected.shape[0])
    return res_annot, selected, summary
