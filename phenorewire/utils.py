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


def save_json(obj: dict, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


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


