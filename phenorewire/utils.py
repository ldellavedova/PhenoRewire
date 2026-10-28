# ── PhenoRewire · Shared utilities ───────────────────────────────────────────
# I/O helpers, directory creation, normalization utilities, and shared functions
# used across modules.
from __future__ import annotations

from pathlib import Path
import json
from typing import Dict, Any, List

import numpy as np
import pandas as pd
import logging
import re

logger = logging.getLogger(__name__)


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def parse_time_numeric(series: pd.Series) -> pd.Series:
    """Times as numbers; "3h", "3 hours" and "3" all parse as 3."""
    s = series.astype(str).str.strip().str.lower()
    s = s.str.replace("hours", "", regex=False).str.replace("hour", "", regex=False).str.replace("h", "", regex=False)
    s = s.str.strip()
    return pd.to_numeric(s, errors="coerce")


def timepoint_label(t: float) -> str:
    """Filename-safe label for a timepoint.

    Integer-valued timepoints keep the plain "3h" form.  Fractional ones keep their
    value with the decimal point spelled out ("0.5" -> "0p5h"), so two distinct
    timepoints can never collapse onto the same directory name.
    """
    t = float(t)
    if t.is_integer():
        return f"{int(t)}h"
    return f"{t:g}".replace("-", "neg").replace(".", "p") + "h"


def safe_log2p1(X: pd.DataFrame) -> pd.DataFrame:
    return np.log2(X.astype(float) + 1.0)


def median_normalize(X: pd.DataFrame) -> pd.DataFrame:
    """Median normalization (column-wise) with safe handling of zero/NaN medians.

    If a column median is 0 (or NaN), its scale factor is set to 1.0 and a warning is logged.
    """
    col_medians = X.median(axis=0, skipna=True)
    target = float(np.nanmedian(col_medians.values))

    scale = target / col_medians.replace(0, np.nan)
    bad = scale.isna()
    if bad.any():
        logger.warning(
            "median_normalize: %d/%d columns have median=0 or NaN; leaving them unscaled (scale=1.0).",
            int(bad.sum()),
            int(len(scale)),
        )
        scale = scale.fillna(1.0)

    return X.mul(scale, axis=1)


PQN_MIN_COMPLETE_FEATURES = 10


def pqn_normalize(X: pd.DataFrame) -> pd.DataFrame:
    """Probabilistic quotient normalization (Dieterle et al., Anal Chem 2006).

    The reference profile is each feature's median over the samples.  A sample's
    dilution factor is the median, over its features, of its intensity divided
    by the reference; the sample is divided by that factor.  Undetected values
    (0) take no part in either median, so sparse samples are normalised on what
    they contain.

    The reference is built from the features detected in every sample when at
    least PQN_MIN_COMPLETE_FEATURES of them exist.  A feature missing from some
    samples has its median taken over a different set of samples than the
    others, and when the missingness follows dilution (a feature undetected in
    the most dilute samples) that biases every quotient.  With too few complete
    features, every feature contributes, its median taken over the samples where
    it was detected.  A sample with no usable quotient is left unscaled and a
    warning is logged.
    """
    detected = X.where(X > 0)
    complete = detected.notna().all(axis=1)
    if int(complete.sum()) >= PQN_MIN_COMPLETE_FEATURES:
        detected = detected.loc[complete]
    reference = detected.median(axis=1, skipna=True)
    quotients = detected.div(reference.where(reference > 0), axis=0)
    factor = quotients.median(axis=0, skipna=True)
    bad = factor.isna() | (factor <= 0)
    if bad.any():
        logger.warning(
            "pqn_normalize: %d/%d samples have no feature detected in common with the reference; "
            "leaving them unscaled.",
            int(bad.sum()),
            int(len(factor)),
        )
        factor = factor.mask(bad, 1.0)
    return X.div(factor, axis=1)


def tic_normalize(X: pd.DataFrame) -> pd.DataFrame:
    """Total-signal normalization: every sample is scaled to the median total intensity."""
    totals = X.sum(axis=0, skipna=True)
    scale = float(np.nanmedian(totals.values)) / totals.replace(0, np.nan)
    bad = scale.isna()
    if bad.any():
        logger.warning("tic_normalize: %d/%d samples have zero total intensity; leaving them unscaled.",
                       int(bad.sum()), int(len(scale)))
        scale = scale.fillna(1.0)
    return X.mul(scale, axis=1)


# How each NORM_METHOD is named in the console and the report.
NORM_LABELS = {
    "median": "median normalisation",
    "pqn": "probabilistic quotient normalisation (PQN)",
    "tic": "total-signal (TIC) normalisation",
    "none": "no normalisation",
}

# NORM_METHOD -> function; "none" keeps the intensities as they are.
NORMALIZERS = {
    "median": median_normalize,
    "pqn": pqn_normalize,
    "tic": tic_normalize,
    "none": lambda X: X,
}


def _json_default(value: Any):
    """Fallback encoder for values pandas/numpy hand back that json cannot serialize.

    Reporting runs last, so an unencodable numpy scalar here would throw away a whole
    pipeline run.  Convert what we can, stringify the rest.
    """
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _json_safe(obj: Any) -> Any:
    """Recursively replace non-finite floats with None.

    ``json.dump`` writes bare ``NaN``/``Infinity`` tokens, which are not valid JSON
    and are rejected by strict parsers (``JSON.parse``, ``jq``).  The reports are
    meant to be machine-readable, so emit ``null`` instead.
    """
    if isinstance(obj, float):
        return obj if np.isfinite(obj) else None
    if isinstance(obj, np.generic):
        return _json_safe(obj.item())
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def save_json(obj: dict, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_json_safe(obj), f, indent=2, default=_json_default)


def assign_groups(meta: pd.DataFrame, group_definition: Dict[str, Dict[str, Any]]) -> List[str]:
    """Assign each sample to the first matching group in GROUP_DEFINITION.

    Robust matching across common CSV/YAML type mismatches:
      - handles int/float/bool vs string mismatches (e.g. "0" vs 0)
      - strips whitespace
      - treats empty strings / NaN as missing
    """

    def _norm(x):
        if x is None:
            return None
        try:
            if pd.isna(x):
                return None
        except Exception:
            pass

        # pandas hands back numpy scalars whenever a row is not object-dtype, and
        # numpy scalars are not instances of the Python int/bool types checked below.
        # Unwrap them first, or a metadata column of plain integers never matches the
        # integer written in the YAML.
        if isinstance(x, np.generic):
            x = x.item()

        if isinstance(x, str):
            s = x.strip()
            if s == "":
                return None
            sl = s.lower()
            if sl in {"true", "t", "yes", "y"}:
                return True
            if sl in {"false", "f", "no", "n"}:
                return False
            if re.fullmatch(r"[+-]?\d+", s):
                try:
                    return int(s)
                except Exception:
                    return s
            if re.fullmatch(r"[+-]?\d+\.\d+", s):
                try:
                    return float(s)
                except Exception:
                    return s
            return s

        if isinstance(x, bool):
            return x
        if isinstance(x, int):
            return x
        if isinstance(x, float):
            return float(x)

        return str(x).strip()

    groups: List[str] = []
    for _, row in meta.iterrows():
        assigned = None
        for group_name, conditions in group_definition.items():
            match = True
            for col, expected in conditions.items():
                if col not in row.index:
                    match = False
                    break
                if _norm(row[col]) != _norm(expected):
                    match = False
                    break
            if match:
                assigned = str(group_name).strip()
                break
        groups.append(assigned or "Unknown")
    return groups


