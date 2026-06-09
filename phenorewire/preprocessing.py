# ── PhenoRewire · Step 1 — Preprocessing ─────────────────────────────────────
# Filters features by global prevalence, applies median normalization across
# samples, and log2(x+1) transforms intensities.
# Input: raw feature table + metadata. Output: filtered intensity matrix per group.
from __future__ import annotations

from pathlib import Path
import re
from typing import Optional, Tuple, Dict, Any

import numpy as np
import pandas as pd
import logging

from .utils import ensure_dir, median_normalize, safe_log2p1, assign_groups

logger = logging.getLogger(__name__)


def prepare_matrices(
    df: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    feature_id_col: str,
    mz_col: str,
    rt_col: str,
    name_col: Optional[str],
    intensity_regex: str,
    meta_sample_col: str,
    outdir: Path,
    min_nonzero_intensity: float,
    min_presence_global: float,
    norm_method: str,
    log_transform: bool,
    group_definition: Dict[str, Dict[str, Any]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ensure_dir(outdir)
    logger.info("Preparing matrices...")

    def extract_sample_name(intensity_col: str) -> str:
        s = str(intensity_col).strip()
        # If INTENSITY_REGEX matches at the end of the column name, strip that suffix.
        m = re.search(intensity_regex, s)
        if m and m.end() == len(s):
            s = s[: m.start()].strip()
        # Common fallbacks (MZmine-style)
        for suffix in ["Peak area", "Peak height", "Area", "Height"]:
            if s.endswith(suffix):
                s = s[: -len(suffix)].strip()
        # Remove trailing separators
        s = re.sub(r"[:\-_\\s]+$", "", s).strip()
        return s

    def norm_sample_id(x: str) -> str:
        x = str(x).strip().replace("\\", "/")
        x = x.split("/")[-1]  # drop path
        # drop common extensions
        for ext in [".mzML", ".mzml", ".mzXML", ".mzxml", ".raw", ".RAW"]:
            if x.endswith(ext):
                x = x[: -len(ext)]
                break
        return x.strip()

    # -------------------------
    # Find intensity columns
    # -------------------------
    intensity_cols = [c for c in df.columns if re.search(intensity_regex, str(c))]
    if not intensity_cols:
        raise ValueError("No intensity columns found. Check INTENSITY_REGEX or your column names.")

    samples_in_matrix_raw = [extract_sample_name(c) for c in intensity_cols]
    samples_in_matrix = [norm_sample_id(s) for s in samples_in_matrix_raw]

    # Build intensity matrix (features x samples)
    X = df[intensity_cols].copy()
    X.columns = samples_in_matrix

    # -------------------------
    # Clean intensities
    # -------------------------
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(float)

    if min_nonzero_intensity > 0:
        X = X.where(X >= float(min_nonzero_intensity), 0.0)

    # presence filter (global)
    presence = (X > 0).mean(axis=1)
    keep = presence >= float(min_presence_global)
    X = X.loc[keep].copy()
    logger.info("Global presence filter kept %d/%d features.", int(keep.sum()), int(len(keep)))

    # -------------------------
    # Normalize + log
    # -------------------------
    if norm_method.lower() == "median":
        X = median_normalize(X)
    else:
        raise ValueError(f"Unsupported norm_method: {norm_method}")

    if log_transform:
        X = safe_log2p1(X)

    # Remove zero-variance features (constant rows produce NaN Spearman r)
    row_std = X.std(axis=1, skipna=True)
    nonzero_var = row_std > 0
    n_dropped = int((~nonzero_var).sum())
    if n_dropped > 0:
        logger.warning(
            "Removed %d zero-variance feature(s) before correlation analysis.",
            n_dropped,
        )
        X = X.loc[nonzero_var].copy()

    # -------------------------
    # Align metadata to matrix samples
    # -------------------------
    meta2 = meta.copy()
    if meta_sample_col not in meta2.columns:
        raise ValueError(f"Metadata is missing META_SAMPLE_COL='{meta_sample_col}'")

    meta2[meta_sample_col] = meta2[meta_sample_col].astype(str).str.strip().map(norm_sample_id)

    # Keep only metadata samples that are present in the matrix (after normalization)
    meta_samples = meta2[meta_sample_col].astype(str).tolist()
    meta_set = set(meta_samples)
    matrix_set = set(samples_in_matrix)

    missing_in_meta = [s for s in samples_in_matrix if s not in meta_set]
    extra_in_meta = [s for s in meta_samples if s not in matrix_set]

    if len(missing_in_meta) == len(samples_in_matrix):
        raise ValueError(
            "No metadata samples match matrix columns (after normalization).\n"
            "Check META_SAMPLE_COL values vs intensity column sample names.\n"
            f"Example matrix samples: {samples_in_matrix[:10]}\n"
            f"Example metadata samples: {meta_samples[:10]}"
        )

    if missing_in_meta:
        raise ValueError(
            "Some matrix samples are missing in metadata (after normalization).\n"
            f"Missing (first 15): {missing_in_meta[:15]}\n"
            "Fix your metadata sample IDs or intensity column naming."
        )

    # Detect duplicates in metadata after normalization
    dup_mask = meta2[meta_sample_col].duplicated(keep=False)
    if dup_mask.any():
        dups = meta2.loc[dup_mask, meta_sample_col].astype(str).unique().tolist()
        raise ValueError(
            "Duplicate sample IDs found in metadata after normalization.\n"
            f"Duplicates (first 15): {dups[:15]}\n"
            "Make sample IDs unique before running the pipeline."
        )

    # Optional warning for extra metadata rows not used
    if extra_in_meta:
        logger.warning(
            "Metadata contains %d samples not present in the matrix (ignored). Example: %s",
            len(set(extra_in_meta)),
            extra_in_meta[:10],
        )

    # Reorder metadata to match matrix sample order
    meta2 = meta2.set_index(meta_sample_col).loc[samples_in_matrix].reset_index()

    # Assign groups
    meta2["group"] = assign_groups(meta2, group_definition)
    meta2.to_csv(outdir / "metadata_aligned.csv", index=False)

    # -------------------------
    # Validate required feature columns
    # -------------------------
    for c in [feature_id_col, mz_col, rt_col]:
        if c not in df.columns:
            raise ValueError(f"Feature table is missing required column '{c}'")

    feat_cols = [feature_id_col, mz_col, rt_col]
    if name_col is not None and name_col in df.columns:
        feat_cols.append(name_col)

    feat_anno = df.loc[X.index, feat_cols].copy()
    feat_anno = feat_anno.rename(
        columns={
            feature_id_col: "feature_id",
            mz_col: "mz",
            rt_col: "rt",
            name_col: "name" if name_col else name_col,
        }
    )

    # Ensure feature_id is string
    feat_anno["feature_id"] = feat_anno["feature_id"].astype(str)

    # Set X index to feature_id (and align to annotation)
    X_feat_by_samp = X.copy()
    X_feat_by_samp.index = feat_anno["feature_id"].astype(str).values
    feat_anno = feat_anno.set_index("feature_id", drop=True)

    # Save diagnostics
    feat_anno.to_csv(outdir / "feature_annotation.csv")
    X_feat_by_samp.to_csv(outdir / "X_corr_matrix.csv")

    return X_feat_by_samp, meta2, feat_anno
