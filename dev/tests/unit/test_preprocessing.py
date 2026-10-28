"""Input guards in prepare_matrices() and the shared helpers it uses.

These cover the parts of the pipeline that read user files, where a malformed
input used to propagate silently into the networks instead of stopping the run.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from phenorewire.preprocessing import prepare_matrices
from phenorewire.utils import assign_groups, median_normalize, save_json

REPO = Path(__file__).resolve().parents[3]

GROUPS = {"Reference": {"group_label": "reference"}, "Case": {"group_label": "case"}}


def _feature_table(n_features: int = 6, n_samples: int = 8) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    data = {
        "feature_id": [f"F{i:03d}" for i in range(n_features)],
        "mz": np.linspace(100.0, 500.0, n_features),
        "rt": np.linspace(1.0, 9.0, n_features),
        "name": [f"compound_{i}" for i in range(n_features)],
    }
    for s in range(n_samples):
        data[f"Sample_{s:02d}"] = rng.uniform(500, 5000, n_features)
    return pd.DataFrame(data)


def _metadata(n_samples: int = 8) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": [f"Sample_{s:02d}" for s in range(n_samples)],
            "group_label": ["reference"] * (n_samples // 2) + ["case"] * (n_samples // 2),
        }
    )


def _run(df: pd.DataFrame, meta: pd.DataFrame, outdir: Path):
    return prepare_matrices(
        df=df,
        meta=meta,
        feature_id_col="feature_id",
        mz_col="mz",
        rt_col="rt",
        name_col="name",
        intensity_regex="^Sample_",
        meta_sample_col="sample_id",
        outdir=outdir,
        min_nonzero_intensity=0.0,
        min_presence_global=0.0,
        norm_method="median",
        log_transform=True,
        group_definition=GROUPS,
    )


def test_happy_path_returns_aligned_matrix(tmp_path: Path) -> None:
    X, meta2, feat_anno = _run(_feature_table(), _metadata(), tmp_path)
    assert X.shape == (6, 8)
    assert list(X.columns) == meta2["sample_id"].tolist()
    assert not X.index.duplicated().any()
    assert set(meta2["group"]) == {"Reference", "Case"}
    assert list(feat_anno.columns) == ["mz", "rt", "name"]


def test_duplicate_feature_ids_are_rejected(tmp_path: Path) -> None:
    """Duplicate IDs used to create self-loops and duplicated .loc[] rows."""
    df = _feature_table()
    df.loc[1, "feature_id"] = df.loc[0, "feature_id"]
    with pytest.raises(ValueError, match="duplicated value"):
        _run(df, _metadata(), tmp_path)


def test_empty_feature_ids_are_rejected(tmp_path: Path) -> None:
    df = _feature_table()
    df.loc[2, "feature_id"] = ""
    with pytest.raises(ValueError, match="empty or missing"):
        _run(df, _metadata(), tmp_path)


def test_sample_missing_from_metadata_is_rejected(tmp_path: Path) -> None:
    meta = _metadata().iloc[:-1]
    with pytest.raises(ValueError, match="missing in metadata"):
        _run(_feature_table(), meta, tmp_path)


def test_duplicate_sample_ids_are_rejected(tmp_path: Path) -> None:
    meta = _metadata()
    meta = pd.concat([meta, meta.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate sample IDs"):
        _run(_feature_table(), meta, tmp_path)


def test_intensity_regex_matching_nothing_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="No intensity columns"):
        prepare_matrices(
            df=_feature_table(),
            meta=_metadata(),
            feature_id_col="feature_id",
            mz_col="mz",
            rt_col="rt",
            name_col="name",
            intensity_regex="^NoSuchPrefix_",
            meta_sample_col="sample_id",
            outdir=tmp_path,
            min_nonzero_intensity=0.0,
            min_presence_global=0.0,
            norm_method="median",
            log_transform=True,
            group_definition=GROUPS,
        )


def test_raw_constant_features_are_dropped(tmp_path: Path) -> None:
    """Median normalization would otherwise turn them into a perfectly correlated clique.

    A row that is constant in the raw data becomes proportional to the per-column
    scale factors after normalization, so a variance check applied only afterwards
    never sees it — and every such row correlates at rho = 1 with every other one.
    """
    df = _feature_table()
    sample_cols = [c for c in df.columns if c.startswith("Sample_")]
    df.loc[3, sample_cols] = 1000.0
    df.loc[4, sample_cols] = 2000.0

    X, _, _ = _run(df, _metadata(), tmp_path)

    assert "F003" not in X.index
    assert "F004" not in X.index
    assert X.shape[0] == 4


def test_dropping_constant_features_leaves_every_other_feature_untouched(tmp_path: Path) -> None:
    """The filter must remove the artefact without shifting anything else.

    Constant features still take part in normalization, so the scale factors stay
    exactly what they were before this filter existed.  Dropping the constant rows
    afterwards therefore cannot change any surviving feature's values.  Dropping
    them *before* normalization would move every column median and shift every
    result, which is a far larger change than removing the artefact requires.
    """
    from phenorewire.utils import median_normalize, safe_log2p1

    sample_cols = [f"Sample_{s:02d}" for s in range(8)]
    df = _feature_table()
    df.loc[3, sample_cols] = 1000.0
    df.loc[4, sample_cols] = 2000.0

    # Reference: the pre-0.1.1 path — normalize with every feature, drop nothing.
    raw = df[sample_cols].astype(float)
    raw.index = df["feature_id"]
    expected = safe_log2p1(median_normalize(raw))

    X, _, _ = _run(df, _metadata(), tmp_path)

    survivors = list(X.index)
    assert survivors == ["F000", "F001", "F002", "F005"]
    # Bit-identical, not almost-equal: the constant rows still set the medians.
    assert np.array_equal(X.to_numpy(), expected.loc[survivors].to_numpy())


def test_constant_features_do_not_create_a_spurious_clique(tmp_path: Path) -> None:
    from phenorewire.networks import correlation_network

    df = _feature_table(n_features=10, n_samples=12)
    sample_cols = [c for c in df.columns if c.startswith("Sample_")]
    for i in range(4):
        df.loc[i, sample_cols] = 1000.0 * (i + 1)

    meta = pd.DataFrame(
        {
            "sample_id": sample_cols,
            "group_label": ["reference"] * 6 + ["case"] * 6,
        }
    )
    X, _, _ = _run(df, meta, tmp_path)
    _, edges, _ = correlation_network(
        X, abs_r_thr=0.6, fdr_alpha=0.05, adaptive_thresholds=False
    )

    dropped = {f"F{i:03d}" for i in range(4)}
    touching_dropped = [
        (a, b) for a, b in zip(edges["node_a"], edges["node_b"]) if a in dropped or b in dropped
    ]
    assert touching_dropped == []


# ---------------------------------------------------------------------------
# assign_groups: CSV/YAML type coercion
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("cell", "expected"),
    [
        (1, 1), ("1", 1), (1.0, 1),
        ("true", True), (True, True),
        ("  case  ", "case"),
    ],
)
def test_assign_groups_matches_across_type_mismatches(cell, expected) -> None:
    meta = pd.DataFrame({"flag": [cell]})
    groups = assign_groups(meta, {"Hit": {"flag": expected}})
    assert groups == ["Hit"]


def test_assign_groups_returns_unknown_when_nothing_matches() -> None:
    meta = pd.DataFrame({"group_label": ["something_else"]})
    assert assign_groups(meta, GROUPS) == ["Unknown"]


def test_assign_groups_takes_the_first_matching_definition() -> None:
    meta = pd.DataFrame({"group_label": ["reference"], "timepoint": [1]})
    definition = {
        "Narrow": {"group_label": "reference", "timepoint": 1},
        "Broad": {"group_label": "reference"},
    }
    assert assign_groups(meta, definition) == ["Narrow"]


def test_assign_groups_requires_every_condition_column_to_exist() -> None:
    meta = pd.DataFrame({"group_label": ["reference"]})
    assert assign_groups(meta, {"G": {"group_label": "reference", "absent": 1}}) == ["Unknown"]


# ---------------------------------------------------------------------------
# Shared utils
# ---------------------------------------------------------------------------

def test_median_normalize_leaves_zero_median_columns_unscaled() -> None:
    X = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [0.0, 0.0, 0.0]})
    out = median_normalize(X)
    assert out["b"].tolist() == [0.0, 0.0, 0.0]
    assert np.isfinite(out.to_numpy()).all()


def test_save_json_survives_numpy_scalars(tmp_path: Path) -> None:
    """Reporting runs last; an unencodable numpy scalar must not discard the run."""
    import json

    payload = {
        "count": np.int64(3),
        "score": np.float64(0.5),
        "flag": np.bool_(True),
        "vector": np.array([1, 2]),
        "missing": np.nan,
    }
    path = tmp_path / "out.json"
    save_json(payload, path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == {
        "count": 3,
        "score": 0.5,
        "flag": True,
        "vector": [1, 2],
        "missing": None,
    }


# ---------------------------------------------------------------------------
# Normalisation methods
# ---------------------------------------------------------------------------

def _diluted(sparse: bool = False) -> tuple[pd.DataFrame, np.ndarray]:
    """One profile, measured in six samples at different dilutions."""
    rng = np.random.default_rng(3)
    profile = rng.lognormal(8, 1, size=40)
    dilution = np.array([1.0, 0.5, 2.0, 0.8, 1.5, 3.0])
    X = pd.DataFrame(np.outer(profile, dilution), columns=[f"S{i}" for i in range(6)])
    if sparse:
        # most features absent in the first three samples
        X.iloc[:25, :3] = 0.0
    return X, dilution


def test_pqn_removes_dilution_even_in_sparse_samples() -> None:
    from phenorewire.utils import median_normalize, pqn_normalize

    for sparse in (False, True):
        X, _ = _diluted(sparse)
        out = pqn_normalize(X)
        detected = X > 0
        ratios = out.where(detected).div(out.where(detected).iloc[:, -1], axis=0).to_numpy()
        # Undetected cells are NaN here. pandas 2 dropped them from stack(), pandas 3
        # keeps them, so select the detected cells explicitly instead of relying on
        # either, and check that none of them went missing.
        finite = np.isfinite(ratios)
        assert int(finite.sum()) == int(detected.to_numpy().sum())
        np.testing.assert_allclose(ratios[finite], 1.0, rtol=1e-12)
    # Where median normalisation cannot: those samples have median 0.
    X, _ = _diluted(sparse=True)
    assert (X.median(axis=0).iloc[:3] == 0).all()
    np.testing.assert_allclose(median_normalize(X).iloc[:, :3].to_numpy(), X.iloc[:, :3].to_numpy())


def test_tic_equalises_totals_and_none_changes_nothing() -> None:
    from phenorewire.utils import NORMALIZERS

    X, _ = _diluted()
    totals = NORMALIZERS["tic"](X).sum(axis=0)
    np.testing.assert_allclose(totals.to_numpy(), totals.iloc[0])
    pd.testing.assert_frame_equal(NORMALIZERS["none"](X), X)


def test_norm_method_none_only_log_transforms(tmp_path: Path) -> None:
    df, meta = _feature_table(), _metadata()
    X, _, _ = prepare_matrices(
        df=df, meta=meta, feature_id_col="feature_id", mz_col="mz", rt_col="rt", name_col="name",
        intensity_regex="^Sample_", meta_sample_col="sample_id", outdir=tmp_path,
        min_nonzero_intensity=0.0, min_presence_global=0.0, norm_method="none", log_transform=True,
        group_definition=GROUPS,
    )
    raw = df.set_index("feature_id").loc[X.index, X.columns].astype(float)
    np.testing.assert_allclose(X.to_numpy(), np.log2(raw.to_numpy() + 1.0))


def test_config_accepts_the_normalisation_methods(tmp_path: Path) -> None:
    from phenorewire.config import LeanConfig

    base = yaml.safe_load((REPO / "config_phenotype.yaml").read_text(encoding="utf-8"))
    for key in ("DATA_MATRIX", "METADATA", "SPECREBOOT_NETWORK"):
        base[key] = str(REPO / base[key])
    for method in ("median", "PQN", "tic", "none"):
        assert LeanConfig(**{**base, "NORM_METHOD": method}).NORM_METHOD == method.lower()
    with pytest.raises(Exception, match="NORM_METHOD must be one of"):
        LeanConfig(**{**base, "NORM_METHOD": "quantile"})
