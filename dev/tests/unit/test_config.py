from pathlib import Path

import pytest

from phenorewire.config import LeanConfig



def _base_config(tmp_path: Path) -> dict:
    data = tmp_path / "data.csv"
    meta = tmp_path / "meta.csv"
    data.write_text("x\n1\n", encoding="utf-8")
    meta.write_text("x\n1\n", encoding="utf-8")
    return {
        "DATA_MATRIX": str(data),
        "METADATA": str(meta),
        "OUTDIR": str(tmp_path / "out"),
        "FEATURE_ID_COL": "feature",
        "MZ_COL": "mz",
        "RT_COL": "rt",
        "INTENSITY_REGEX": "Peak area$",
        "META_SAMPLE_COL": "sample",
        "GROUP_DEFINITION": {"A": {"group": 1}, "B": {"group": 0}},
    }


def test_temporal_mode_requires_enabled_temporal_block(tmp_path: Path) -> None:
    cfg = _base_config(tmp_path)
    cfg["ANALYSIS_MODE"] = "temporal"
    with pytest.raises(ValueError, match="requests temporal analysis"):
        LeanConfig(**cfg)


def test_temporal_only_accepts_single_group(tmp_path: Path) -> None:
    cfg = _base_config(tmp_path)
    cfg["ANALYSIS_MODE"] = "temporal"
    cfg["GROUP_DEFINITION"] = {"Copresence": {"group": 1}}
    cfg["TEMPORAL_CORRELATION"] = {
        "enabled": True,
        "phenotype": "Copresence",
        "time_column": "timepoint",
        "n_permutations": 10,
        "fdr_alpha": 0.05,
        "min_samples": 3,
    }
    model = LeanConfig(**cfg)
    assert model.ANALYSIS_MODE == "temporal"
    assert list(model.GROUP_DEFINITION.keys()) == ["Copresence"]


def test_adaptive_selection_defaults_to_false(tmp_path: Path) -> None:
    model = LeanConfig(**_base_config(tmp_path))
    assert model.ADAPTIVE_SELECTION_THRESHOLDS is False


def test_new_config_params_have_correct_defaults(tmp_path: Path) -> None:
    model = LeanConfig(**_base_config(tmp_path))
    assert model.MIN_FEATURES_HARD_STOP == 5
    assert model.TRIAGE_TOPOLOGY_WEIGHT == 0.7
    assert model.TRIAGE_SELECTION_WEIGHT == 0.3
    assert model.SIGN_SWITCH_MIN_R == 0.5
    assert model.LOUVAIN_STABILITY_CHECK is False
    assert model.NETWORK_MIN_EDGES_HARD_STOP == 10
    assert model.NETWORK_MIN_NODES_HARD_STOP == 5


def test_triage_weights_must_sum_to_one(tmp_path: Path) -> None:
    cfg = _base_config(tmp_path)
    cfg["TRIAGE_TOPOLOGY_WEIGHT"] = 0.6
    cfg["TRIAGE_SELECTION_WEIGHT"] = 0.6   # sum = 1.2
    with pytest.raises(ValueError, match="must equal 1.0"):
        LeanConfig(**cfg)


def test_triage_weights_custom_valid_combination(tmp_path: Path) -> None:
    cfg = _base_config(tmp_path)
    cfg["TRIAGE_TOPOLOGY_WEIGHT"] = 0.5
    cfg["TRIAGE_SELECTION_WEIGHT"] = 0.5
    model = LeanConfig(**cfg)
    assert model.TRIAGE_TOPOLOGY_WEIGHT == 0.5
    assert model.TRIAGE_SELECTION_WEIGHT == 0.5


def test_sign_switch_min_r_must_be_in_unit_interval(tmp_path: Path) -> None:
    cfg = _base_config(tmp_path)
    cfg["SIGN_SWITCH_MIN_R"] = 1.5
    with pytest.raises(ValueError):
        LeanConfig(**cfg)


def test_min_features_hard_stop_must_be_positive(tmp_path: Path) -> None:
    cfg = _base_config(tmp_path)
    cfg["MIN_FEATURES_HARD_STOP"] = 0
    with pytest.raises(ValueError):
        LeanConfig(**cfg)
