"""The input check: what it reports, and what it refuses."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenorewire.checks import InputDataError, check_inputs
from phenorewire.config import LeanConfig
from phenorewire.run import run

REPO = Path(__file__).resolve().parents[3]


def _tables(n_per_group=15, n_feat=12, seed=0):
    rng = np.random.default_rng(seed)
    samples = [f"S{i:02d}" for i in range(2 * n_per_group)]
    X = np.exp(rng.normal(8, 1, size=(n_feat, len(samples))))
    feat = pd.DataFrame({"feature_id": [f"F{i}" for i in range(n_feat)], "mz": 100.0, "rt": 1.0})
    feat = pd.concat([feat, pd.DataFrame(X, columns=samples)], axis=1)
    meta = pd.DataFrame({"sample": samples, "group": ["A"] * n_per_group + ["B"] * n_per_group,
                         "time": [0, 4] * n_per_group})
    return feat, meta


def _config(tmp_path, feat, meta, **extra):
    feat.to_csv(tmp_path / "f.csv", index=False)
    meta.to_csv(tmp_path / "m.csv", index=False)
    return LeanConfig(**{
        "DATA_MATRIX": str(tmp_path / "f.csv"), "METADATA": str(tmp_path / "m.csv"),
        "OUTDIR": str(tmp_path / "out"), "ANALYSIS_MODE": "phenotype",
        "GROUP_DEFINITION": {"A": {"group": "A"}, "B": {"group": "B"}},
        "FEATURE_ID_COL": "feature_id", "MZ_COL": "mz", "RT_COL": "rt",
        "INTENSITY_REGEX": "^S", "META_SAMPLE_COL": "sample",
        **extra,
    })


def test_clean_data_pass_and_describe_the_design(tmp_path):
    feat, meta = _tables()
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert chk.ok, chk.errors
    assert not chk.warnings
    text = chk.text()
    assert "12 features x 30 samples" in text
    assert "A (15), B (15)" in text
    assert "within the operating range" in text
    assert chk.comparisons[0]["n"] == [15, 15]


def test_bundled_example_passes():
    chk = check_inputs(LeanConfig.from_yaml(REPO / "config_phenotype.yaml"))
    assert chk.ok and not chk.warnings


def test_autoscaled_data_are_refused(tmp_path):
    """Mean-centred, unit-variance values: half of them negative."""
    feat, meta = _tables()
    cols = [c for c in feat.columns if c.startswith("S")]
    X = feat[cols].to_numpy()
    feat[cols] = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert not chk.ok
    assert "negative" in chk.errors[0] and "autoscaled" in chk.errors[0]


def test_run_stops_on_autoscaled_data_before_any_analysis(tmp_path):
    feat, meta = _tables()
    cols = [c for c in feat.columns if c.startswith("S")]
    feat[cols] = feat[cols] - feat[cols].to_numpy().mean()
    cfg = _config(tmp_path, feat, meta)
    with pytest.raises(InputDataError, match="negative"):
        run(cfg)
    out = Path(cfg.OUTDIR)
    assert (out / "run_info" / "data_check.txt").exists()
    assert not (out / "A_vs_B").exists()


def test_boolean_features_are_named(tmp_path):
    """Excel booleans exported in another language arrive as text."""
    feat, meta = _tables()
    cols = [c for c in feat.columns if c.startswith("S")]
    feat[cols] = feat[cols].astype(object)
    feat.loc[0, cols] = ["VRAI", "FAUX"] * 15
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert chk.ok
    assert any("1 feature(s) hold text" in w and "'VRAI'" in w for w in chk.warnings)


def test_log_scale_data_are_flagged(tmp_path):
    feat, meta = _tables()
    cols = [c for c in feat.columns if c.startswith("S")]
    feat[cols] = np.log2(feat[cols])
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert chk.ok
    assert any("LOG_TRANSFORM: false" in w for w in chk.warnings)
    assert not check_inputs(_config(tmp_path, feat, meta, LOG_TRANSFORM=False)).warnings


def test_samples_missing_from_metadata_are_listed(tmp_path):
    feat, meta = _tables()
    chk = check_inputs(_config(tmp_path, feat, meta.iloc[2:]))
    assert not chk.ok
    assert "S00, S01" in chk.errors[0]


def test_unmatched_sample_names_show_both_sides(tmp_path):
    feat, meta = _tables()
    meta["sample"] = "X_" + meta["sample"]
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert "No sample of the feature table has a row in the metadata" in chk.errors[0]
    assert "S00" in chk.errors[0] and "X_S00" in chk.errors[0]


def test_wrong_group_column_is_named(tmp_path):
    feat, meta = _tables()
    cfg = _config(tmp_path, feat, meta, GROUP_DEFINITION={"A": {"Group": "A"}, "B": {"Group": "B"}})
    chk = check_inputs(cfg)
    assert "do not exist: Group." in chk.errors[0]
    assert "sample, group, time" in chk.errors[0]


def test_samples_matching_no_group_are_named(tmp_path):
    """QC injections and blanks, and any mistyped group value, must not drop silently.

    assign_groups labels an unmatched sample "Unknown" rather than leaving it null,
    so a nullness test here would never fire.
    """
    feat, meta = _tables()
    for name in ("S90", "S91"):
        feat[name] = feat["S00"]
        meta = pd.concat([meta, pd.DataFrame([{"sample": name, "group": "QC", "time": 0}])],
                         ignore_index=True)
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert chk.ok
    warning = next(w for w in chk.warnings if "match no group" in w)
    assert "2 sample(s)" in warning
    assert "S90" in warning and "S91" in warning


def test_small_groups_warn_with_operating_range(tmp_path):
    feat, meta = _tables(n_per_group=8)
    chk = check_inputs(_config(tmp_path, feat, meta))
    assert chk.ok
    assert any("below the operating range" in w for w in chk.warnings)


def test_temporal_plan_and_missing_timepoint(tmp_path):
    feat, meta = _tables()
    temporal = dict(
        ANALYSIS_MODE="temporal",
        TEMPORAL_CORRELATION={"enabled": True, "phenotype": "A", "time_column": "time"},
        TIMEPOINT_T1=0, TIMEPOINT_T2=4,
    )
    chk = check_inputs(_config(tmp_path, feat, meta, **temporal))
    assert chk.ok, chk.errors
    assert "timepoints in A: 0h (8), 4h (7)" in chk.text()
    assert chk.comparisons[0]["kind"] == "temporal"

    chk = check_inputs(_config(tmp_path, feat, meta, **{**temporal, "TIMEPOINT_T2": 10}))
    assert not chk.ok
    assert "TIMEPOINT_T2 = 10 has no samples in A; available timepoints: 0h, 4h" in chk.errors[0]


def test_sparse_samples_suggest_pqn(tmp_path):
    feat, meta = _tables(n_feat=12)
    cols = [c for c in feat.columns if c.startswith("S")]
    feat.loc[:7, cols[:3]] = 0.0          # 8 of 12 features absent in three samples
    chk = check_inputs(_config(tmp_path, feat, meta, MIN_PRESENCE_GLOBAL=0.0))
    assert any("NORM_METHOD: pqn" in w and "3 sample(s)" in w for w in chk.warnings)
    assert "normalisation: median normalisation (NORM_METHOD: median), then log2(x + 1)" in chk.text()
    pqn = check_inputs(_config(tmp_path, feat, meta, MIN_PRESENCE_GLOBAL=0.0, NORM_METHOD="pqn"))
    assert not any("NORM_METHOD: pqn" in w for w in pqn.warnings)
