# ── PhenoRewire · Configuration & validation ─────────────────────────────────
# Validates the YAML config file using Pydantic v2.
# Called first by run.py before any data is loaded.
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Any
import re

import yaml
from pydantic import BaseModel, Field, ConfigDict
from pydantic import field_validator, model_validator


# ============================================================
# Sub-configs
# ============================================================

class TemporalCorrelationConfig(BaseModel):
    """
    Temporal feature selection within a single phenotype group.
    """
    model_config = ConfigDict(extra="forbid")

    enabled: bool = False
    phenotype: Optional[str] = None
    time_column: Optional[str] = None
    n_permutations: int = 500
    fdr_alpha: float = 0.05
    min_samples: int = 3

    @field_validator("n_permutations")
    @classmethod
    def _v_nperm(cls, v: int) -> int:
        v = int(v)
        if v < 10:
            raise ValueError("TemporalCorrelationConfig.n_permutations must be >= 10.")
        return v

    @field_validator("fdr_alpha")
    @classmethod
    def _v_alpha(cls, v: float) -> float:
        v = float(v)
        if not (0 < v <= 1):
            raise ValueError("TemporalCorrelationConfig.fdr_alpha must be in (0, 1].")
        return v

    @field_validator("min_samples")
    @classmethod
    def _v_min_samples(cls, v: int) -> int:
        v = int(v)
        if v < 3:
            raise ValueError("TemporalCorrelationConfig.min_samples must be >= 3.")
        return v

    @field_validator("phenotype", "time_column")
    @classmethod
    def _v_strip_optional(cls, v):
        if v is None:
            return None
        vv = str(v).strip()
        return vv if vv else None


# ============================================================
# Main Config
# ============================================================

class LeanConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    # Inputs/Outputs
    DATA_MATRIX: Path
    METADATA: Path
    OUTDIR: Path
    ANALYSIS_MODE: str = "auto"

    # Groups
    GROUP_DEFINITION: Dict[str, Dict[str, Any]]

    # Force phenotype direction (A): REF=0, CASE=1
    PHENO_GROUP_REF: Optional[str] = None
    PHENO_GROUP_CASE: Optional[str] = None

    # Temporal rewiring (B): explicit timepoints to compare (e.g. 3 vs 6)
    TIMEPOINT_T1: float = 3
    TIMEPOINT_T2: float = 6

    # Feature table columns
    FEATURE_ID_COL: str
    MZ_COL: str
    RT_COL: str
    NAME_COL: Optional[str] = None
    INTENSITY_REGEX: str

    # Metadata columns
    META_SAMPLE_COL: str
    META_TIME_COL: Optional[str] = None  # fallback for TEMPORAL_CORRELATION.time_column

    # Preprocessing
    MIN_NONZERO_INTENSITY: float = 0.0
    MIN_PRESENCE_GLOBAL: float = 0.25
    LOG_TRANSFORM: bool = True
    NORM_METHOD: str = "median"

    # Phenotype selection
    FDR_ALPHA: float = 0.05
    N_PERMUTATIONS_PHENO: int = 500
    # ADAPTIVE_SELECTION_THRESHOLDS: opt-in relaxation — false by default so reviewers see the
    # base-alpha result; set true only when you have a small cohort and accept the caveat.
    ADAPTIVE_SELECTION_THRESHOLDS: bool = False
    SELECTION_MIN_FEATURES: int = 5
    SELECTION_FDR_ALPHA_CEILING: float = 0.1
    # Hard stop: pipeline raises an error if fewer than this many features pass selection
    # (even after adaptive relaxation).  Prevents silently producing triage from 1-2 features.
    MIN_FEATURES_HARD_STOP: int = 5

    # Correlation network
    CORR_ABS_THRESHOLD: float = 0.6
    CORR_FDR_ALPHA: float = 0.05
    REPORT_TOP_N: int = 15
    ADAPTIVE_NETWORK_THRESHOLDS: bool = True
    NETWORK_MIN_EDGES: int = 3
    NETWORK_MAX_DENSITY: float = 0.25
    NETWORK_ABS_THRESHOLD_FLOOR: float = 0.5
    NETWORK_FDR_ALPHA_CEILING: float = 0.2
    # Hard stops: pipeline raises an error when the final network is too sparse to be meaningful.
    NETWORK_MIN_EDGES_HARD_STOP: int = 10
    NETWORK_MIN_NODES_HARD_STOP: int = 5

    # Sign-switch quality filter
    # Both edges in a shared pair must have |r| >= this threshold to be flagged as a sign switch.
    # Set to 0.0 to disable the filter and flag all sign changes regardless of effect size.
    SIGN_SWITCH_MIN_R: float = 0.5

    # Community detection stability (opt-in)
    # When true: Louvain is run 10 times with different seeds and adjusted mutual information
    # (AMI) is reported.  Low AMI (<0.8) triggers a warning in the final report.
    LOUVAIN_STABILITY_CHECK: bool = False

    # Triage score weights — must sum to 1.0
    # topology_priority_score = mean percentile rank of [degree, betweenness, eigenvector, weighted_degree]
    # selection_support_score = 1 - q_value (strength of phenotype/time association)
    TRIAGE_TOPOLOGY_WEIGHT: float = 0.7
    TRIAGE_SELECTION_WEIGHT: float = 0.3

    # Reproducibility
    RANDOM_SEED: int = 42

    # 2D network annotation expansion (used with --run-2dnetwork)
    SPECREBOOT_NETWORK: Optional[str] = None
    PHENOREWIRE_NETWORK: Optional[str] = None
    MN_HOP_DEPTH: int = 2
    MN_MIN_COSINE: float = 0.70
    ANNOTATION_TOP_N: int = 20

    # Modules (NO mutable defaults!)
    TEMPORAL_CORRELATION: TemporalCorrelationConfig = Field(default_factory=TemporalCorrelationConfig)

    # ------------------------------
    # Field validators
    # ------------------------------

    @field_validator("DATA_MATRIX", "METADATA", "OUTDIR", mode="before")
    @classmethod
    def _coerce_path(cls, v):
        return Path(v) if not isinstance(v, Path) else v

    @field_validator("DATA_MATRIX", "METADATA")
    @classmethod
    def _ensure_inputs_exist(cls, v: Path) -> Path:
        v = Path(v)
        if not v.exists():
            raise FileNotFoundError(f"Input path does not exist: {v}")
        return v

    @field_validator("OUTDIR")
    @classmethod
    def _normalize_outdir(cls, v: Path) -> Path:
        # Validation must not touch the filesystem: a config that fails on a later
        # field would otherwise leave an empty output directory behind.  run()
        # creates OUTDIR when it actually starts.
        v = Path(v)
        if v.exists() and not v.is_dir():
            raise ValueError(f"OUTDIR exists but is not a directory: {v}")
        return v

    @field_validator("ANALYSIS_MODE")
    @classmethod
    def _v_analysis_mode(cls, v: str) -> str:
        vv = str(v).strip().lower()
        allowed = {"auto", "phenotype", "temporal", "both"}
        if vv not in allowed:
            raise ValueError(f"ANALYSIS_MODE must be one of {sorted(allowed)}.")
        return vv

    @field_validator("NORM_METHOD")
    @classmethod
    def _check_norm(cls, v: str) -> str:
        vv = str(v).strip().lower()
        if vv != "median":
            raise ValueError("Only NORM_METHOD='median' is supported.")
        return vv

    @field_validator("MIN_PRESENCE_GLOBAL")
    @classmethod
    def _v_presence(cls, v: float) -> float:
        v = float(v)
        if not (0 <= v <= 1):
            raise ValueError("MIN_PRESENCE_GLOBAL must be in [0, 1].")
        return v

    @field_validator("FDR_ALPHA", "CORR_FDR_ALPHA", "SELECTION_FDR_ALPHA_CEILING")
    @classmethod
    def _v_alpha_main(cls, v: float) -> float:
        v = float(v)
        if not (0 < v <= 1):
            raise ValueError("FDR alpha must be in (0, 1].")
        return v

    @field_validator("CORR_ABS_THRESHOLD")
    @classmethod
    def _v_corr_thr(cls, v: float) -> float:
        v = float(v)
        if not (0 <= v <= 1):
            raise ValueError("CORR_ABS_THRESHOLD must be in [0, 1].")
        return v

    @field_validator("N_PERMUTATIONS_PHENO")
    @classmethod
    def _v_nperm_pheno(cls, v: int) -> int:
        v = int(v)
        if v < 10:
            raise ValueError("N_PERMUTATIONS_PHENO must be >= 10 for stability.")
        return v

    @field_validator(
        "REPORT_TOP_N", "NETWORK_MIN_EDGES", "SELECTION_MIN_FEATURES",
        "MIN_FEATURES_HARD_STOP", "NETWORK_MIN_EDGES_HARD_STOP", "NETWORK_MIN_NODES_HARD_STOP",
    )
    @classmethod
    def _v_positive_int(cls, v: int) -> int:
        v = int(v)
        if v < 1:
            raise ValueError("Value must be >= 1.")
        return v

    @field_validator("RANDOM_SEED")
    @classmethod
    def _v_seed(cls, v: int) -> int:
        # 0 is a perfectly ordinary seed; only negative values are rejected.
        v = int(v)
        if v < 0:
            raise ValueError("RANDOM_SEED must be >= 0.")
        return v

    @field_validator("NETWORK_MAX_DENSITY", "NETWORK_ABS_THRESHOLD_FLOOR", "NETWORK_FDR_ALPHA_CEILING")
    @classmethod
    def _v_unit_interval(cls, v: float) -> float:
        v = float(v)
        if not (0 < v <= 1):
            raise ValueError("Adaptive network parameters must be in (0, 1].")
        return v

    @field_validator("SIGN_SWITCH_MIN_R")
    @classmethod
    def _v_sign_switch_min_r(cls, v: float) -> float:
        v = float(v)
        if not (0.0 <= v < 1.0):
            raise ValueError("SIGN_SWITCH_MIN_R must be in [0, 1).")
        return v

    @field_validator("TRIAGE_TOPOLOGY_WEIGHT", "TRIAGE_SELECTION_WEIGHT")
    @classmethod
    def _v_triage_weight(cls, v: float) -> float:
        v = float(v)
        if not (0.0 <= v <= 1.0):
            raise ValueError("Triage weights must be in [0, 1].")
        return v

    @field_validator("FEATURE_ID_COL", "MZ_COL", "RT_COL", "INTENSITY_REGEX", "META_SAMPLE_COL")
    @classmethod
    def _v_required_str(cls, v: str) -> str:
        vv = str(v).strip()
        if not vv:
            raise ValueError("Required string field cannot be empty.")
        return vv

    @field_validator("INTENSITY_REGEX")
    @classmethod
    def _v_intensity_regex_compiles(cls, v: str) -> str:
        try:
            re.compile(v)
        except re.error as e:
            raise ValueError(f"INTENSITY_REGEX is not a valid regex: {e}")
        return v

    @field_validator("NAME_COL", "META_TIME_COL", "PHENO_GROUP_REF", "PHENO_GROUP_CASE", mode="before")
    @classmethod
    def _v_optional_str(cls, v):
        if v is None:
            return None
        vv = str(v).strip()
        return vv if vv else None

    # ------------------------------
    # Model validators (cross-field)
    # ------------------------------

    @model_validator(mode="after")
    def _check_groups_and_direction(self) -> "LeanConfig":
        # GROUP_DEFINITION
        if not isinstance(self.GROUP_DEFINITION, dict) or len(self.GROUP_DEFINITION) < 1:
            raise ValueError("GROUP_DEFINITION must contain at least one group.")

        cleaned: Dict[str, Dict[str, Any]] = {}
        for k, v in self.GROUP_DEFINITION.items():
            kk = str(k).strip()
            if not kk:
                raise ValueError("GROUP_DEFINITION contains an empty group name.")
            if not isinstance(v, dict) or len(v) == 0:
                raise ValueError(f"GROUP_DEFINITION['{kk}'] must be a non-empty dict of conditions.")
            cleaned[kk] = v

        object.__setattr__(self, "GROUP_DEFINITION", cleaned)

        # Direction
        ref = self.PHENO_GROUP_REF
        case = self.PHENO_GROUP_CASE

        if (ref is None) ^ (case is None):
            raise ValueError("Set both PHENO_GROUP_REF and PHENO_GROUP_CASE, or neither.")

        if ref is not None and case is not None:
            if ref == case:
                raise ValueError("PHENO_GROUP_REF and PHENO_GROUP_CASE must be different.")
            keys = set(cleaned.keys())
            if ref not in keys or case not in keys:
                raise ValueError(
                    "PHENO_GROUP_REF/CASE must be keys in GROUP_DEFINITION. "
                    f"Got ref='{ref}', case='{case}', keys={sorted(keys)}"
                )

        # timepoints
        try:
            t1 = float(self.TIMEPOINT_T1)
            t2 = float(self.TIMEPOINT_T2)
        except Exception:
            raise ValueError("TIMEPOINT_T1 and TIMEPOINT_T2 must be numeric (e.g., 3 and 6).")

        if t1 == t2:
            raise ValueError(
                f"TIMEPOINT_T1 and TIMEPOINT_T2 must be different (both are {t1})."
            )
        if t1 > t2:
            raise ValueError(
                f"TIMEPOINT_T1 ({t1}) must be less than TIMEPOINT_T2 ({t2}). "
                "Swap the values or correct your config."
            )

        object.__setattr__(self, "TIMEPOINT_T1", t1)
        object.__setattr__(self, "TIMEPOINT_T2", t2)

        return self

    @model_validator(mode="after")
    def _check_temporal_config(self) -> "LeanConfig":
        tc = self.TEMPORAL_CORRELATION
        mode = self.ANALYSIS_MODE

        # ergonomic fallback: allow META_TIME_COL to populate tc.time_column if missing
        if tc.enabled and (not tc.time_column) and self.META_TIME_COL:
            object.__setattr__(tc, "time_column", self.META_TIME_COL)

        if tc.enabled:
            if not tc.phenotype or not tc.time_column:
                raise ValueError("TEMPORAL_CORRELATION.enabled=True requires phenotype and time_column.")
            if tc.n_permutations < 10:
                raise ValueError("TEMPORAL_CORRELATION.n_permutations must be >= 10 for stability.")
            if tc.min_samples < 3:
                raise ValueError("TEMPORAL_CORRELATION.min_samples must be >= 3.")
        elif mode in {"temporal", "both"}:
            raise ValueError(
                "ANALYSIS_MODE requests temporal analysis, but TEMPORAL_CORRELATION.enabled is False."
            )

        if mode in {"phenotype", "both"} and len(self.GROUP_DEFINITION) < 2:
            raise ValueError(
                "Phenotype analysis requires at least two groups in GROUP_DEFINITION."
            )
        return self

    @model_validator(mode="after")
    def _check_triage_weights(self) -> "LeanConfig":
        total = self.TRIAGE_TOPOLOGY_WEIGHT + self.TRIAGE_SELECTION_WEIGHT
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                f"TRIAGE_TOPOLOGY_WEIGHT + TRIAGE_SELECTION_WEIGHT must equal 1.0 "
                f"(got {self.TRIAGE_TOPOLOGY_WEIGHT:.4f} + {self.TRIAGE_SELECTION_WEIGHT:.4f} = {total:.6f})."
            )
        return self

    # ------------------------------
    # Loader
    # ------------------------------
    @classmethod
    def from_yaml(cls, path: Path | str) -> "LeanConfig":
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config YAML not found: {path}")

        base_dir = path.parent.resolve()

        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}

        if not isinstance(data, dict):
            raise ValueError("Top-level YAML must be a mapping (key: value).")

        # Resolve relative paths with respect to YAML location (reproducible)
        for key in ("DATA_MATRIX", "METADATA", "OUTDIR", "SPECREBOOT_NETWORK", "PHENOREWIRE_NETWORK"):
            if key in data and data[key] is not None:
                p = Path(str(data[key]))
                if not p.is_absolute():
                    data[key] = str((base_dir / p).resolve())

        return cls(**data)
