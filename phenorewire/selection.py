# ── PhenoRewire · Step 2 — Shared feature-selection core ─────────────────────
# The phenotype and temporal selections differ only in the target they correlate
# against, the names of their output columns and files, and the extra effect-size
# columns the phenotype mode reports.  Everything else — permutation testing, BH
# correction, adaptive relaxation, the hard stop, annotation attachment — lives here
# so the two modes cannot drift apart.
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import logging

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .stats import adaptive_fdr_selection, spearman_permutation_test
from .utils import save_json

logger = logging.getLogger(__name__)

# Columns that may carry a human-readable compound name, in priority order.
_NAME_ALIASES = ("NAME", "compound_name", "metabolite_name", "consensus_annotation", "annotation")


def attach_annotation(res: pd.DataFrame, feat_anno: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Join feature annotation onto a result table and guarantee a usable 'name'.

    ``feat_anno`` may be None, empty, indexed by feature_id, or carry feature_id as
    a column under any of several aliases.  The returned frame always has a 'name'
    column, falling back to the feature_id when no real name is available.
    """
    res = res.copy()
    res.index = res.index.astype(str).str.strip()

    if feat_anno is None or feat_anno.empty:
        anno = pd.DataFrame(index=res.index)
        anno.index.name = "feature_id"
    else:
        anno = feat_anno.copy()
        if anno.index.name != "feature_id":
            if "feature_id" in anno.columns:
                anno["feature_id"] = anno["feature_id"].astype(str).str.strip()
                anno = anno.set_index("feature_id", drop=False)
            else:
                anno.index = anno.index.astype(str).str.strip()
                anno.index.name = "feature_id"

        if "name" not in anno.columns:
            for alt in _NAME_ALIASES:
                if alt in anno.columns:
                    anno = anno.rename(columns={alt: "name"})
                    break

    out = res.join(anno, how="left")

    fallback = pd.Series(out.index.astype(str), index=out.index, name="name")
    if "name" not in out.columns:
        out["name"] = fallback
    else:
        cleaned = out["name"].replace({None: np.nan})
        cleaned = cleaned.where(cleaned.notna(), fallback)
        # Avoid the literal string "nan" leaking in from a stringified missing value.
        cleaned = cleaned.astype(str).replace({"nan": np.nan, "": np.nan})
        out["name"] = cleaned.fillna(fallback)

    return out


def _hard_stop_message(
    *,
    mode: str,
    n_selected: int,
    applied_alpha: float,
    min_features_hard_stop: int,
    adaptive_fdr: bool,
    fdr_alpha_ceiling: float,
    remedy: str,
) -> str:
    adaptive_note = (
        f"  Adaptive FDR relaxation was enabled; ceiling alpha = {fdr_alpha_ceiling}.\n"
        if adaptive_fdr
        else "  ADAPTIVE_SELECTION_THRESHOLDS is False; enable it or increase sample size.\n"
    )
    return (
        f"{mode} feature selection hard stop: only {n_selected} feature(s) passed FDR "
        f"threshold (q <= {applied_alpha:.4f}), fewer than "
        f"MIN_FEATURES_HARD_STOP = {min_features_hard_stop}.\n"
        + adaptive_note
        + f"  Consider: {remedy}, a higher FDR alpha, or enabling "
        "ADAPTIVE_SELECTION_THRESHOLDS."
    )


def run_permutation_selection(
    X: pd.DataFrame,
    target: np.ndarray,
    feat_anno: Optional[pd.DataFrame],
    outdir: Path,
    *,
    mode: str,
    r_col: str,
    q_col: str,
    all_filename: str,
    selected_filename: str,
    summary_filename: str,
    n_permutations: int,
    fdr_alpha: float,
    adaptive_fdr: bool,
    min_selected_features: int,
    fdr_alpha_ceiling: float,
    seed: int,
    min_features_hard_stop: int,
    remedy: str,
    extra_columns: Optional[dict] = None,
    extra_summary: Optional[dict] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Correlate every feature against ``target``, FDR-correct, and select.

    Shared by the phenotype and temporal modes; see those wrappers for the
    mode-specific contract.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    feats = X.index.astype(str).str.strip().tolist()
    mat = X.values.astype(float)

    n_perm = int(n_permutations)
    if n_perm < 10:
        raise ValueError("n_permutations must be >= 10 for stability.")

    r_obs, p_emp = spearman_permutation_test(
        mat, np.asarray(target, dtype=float), n_permutations=n_perm, seed=int(seed)
    )

    # With method="fdr_bh" the returned q-values do not depend on alpha; it only
    # sets the (unused) reject array.  Selection is done against fdr_alpha below.
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
        raise ValueError(
            _hard_stop_message(
                mode=mode,
                n_selected=n_selected,
                applied_alpha=float(selection_meta["applied_fdr_alpha"]),
                min_features_hard_stop=int(min_features_hard_stop),
                adaptive_fdr=bool(adaptive_fdr),
                fdr_alpha_ceiling=float(fdr_alpha_ceiling),
                remedy=remedy,
            )
        )

    columns = {
        "feature_id": feats,
        r_col: r_obs.astype(float),
        "p_empirical": p_emp.astype(float),
        q_col: qvals.astype(float),
        "significant": rej.astype(bool),
    }
    columns.update(extra_columns or {})
    res = pd.DataFrame(columns).set_index("feature_id")

    res_annot = attach_annotation(res, feat_anno)
    res_annot.to_csv(outdir / all_filename, index=True)

    selected = res_annot[res_annot["significant"]].copy()
    selected.to_csv(outdir / selected_filename, index=True)

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
    summary.update(extra_summary or {})
    save_json(summary, outdir / summary_filename)

    logger.info("%s selection done. Selected %d features.", mode, selected.shape[0])
    return res_annot, selected, summary
