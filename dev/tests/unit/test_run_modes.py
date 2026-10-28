"""Mode resolution in run(): an analysis the user asked for must never be skipped silently.

Regression tests for the case where ANALYSIS_MODE requests a phenotype contrast
but the data does not yield exactly two groups.  Before, run() logged an INFO
line, produced no analysis and exited 0 with a report claiming success.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from phenorewire.config import LeanConfig
from phenorewire.run import run

REPO_ROOT = Path(__file__).resolve().parents[3]
BASE_CONFIG = REPO_ROOT / "config_phenotype.yaml"

# Three groups, all populated: reference splits by timepoint, case does not.
THREE_GROUPS = {
    "RefT1": {"group_label": "reference", "timepoint": 1},
    "RefT2": {"group_label": "reference", "timepoint": 2},
    "Case": {"group_label": "case"},
}


def _write_config(tmp_path: Path, *, n_permutations: int = 20, **overrides) -> Path:
    """Write a config derived from the bundled phenotype template.

    ``n_permutations`` defaults low because most tests here assert that run()
    raises before selection ever happens.  Tests that need selection to succeed
    must pass enough permutations for the smallest achievable p-value,
    1/(n+1), to survive BH correction.
    """
    data = copy.deepcopy(yaml.safe_load(BASE_CONFIG.read_text(encoding="utf-8")))
    # Resolve the bundled data paths against the repo root, not tmp_path.
    for key in ("DATA_MATRIX", "METADATA"):
        data[key] = str((REPO_ROOT / data[key]).resolve())
    data["OUTDIR"] = str(tmp_path / "out")
    data["N_PERMUTATIONS_PHENO"] = int(n_permutations)
    data.update(overrides)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


def test_phenotype_mode_raises_when_groups_are_not_a_clean_pair(tmp_path: Path) -> None:
    cfg = _write_config(
        tmp_path,
        ANALYSIS_MODE="phenotype",
        GROUP_DEFINITION=THREE_GROUPS,
        PHENO_GROUP_REF=None,
        PHENO_GROUP_CASE=None,
    )
    with pytest.raises(ValueError, match="two groups"):
        run(LeanConfig.from_yaml(cfg))


def test_explicit_ref_case_restricts_the_contrast_instead_of_refusing(tmp_path: Path) -> None:
    """Naming the two groups should run that contrast, not abort on the third group."""
    cfg = _write_config(
        tmp_path,
        n_permutations=1000,
        ANALYSIS_MODE="phenotype",
        GROUP_DEFINITION=THREE_GROUPS,
        PHENO_GROUP_REF="RefT1",
        PHENO_GROUP_CASE="Case",
    )
    config = LeanConfig.from_yaml(cfg)
    run(config)

    outdir = Path(config.OUTDIR)
    report_json = json.loads((outdir / "report.json").read_text(encoding="utf-8"))
    assert report_json["phenotype"]["ref_group"] == "RefT1"
    assert report_json["phenotype"]["case_group"] == "Case"

    # Only the two contrasted groups may contribute samples.
    used = pd.read_csv(outdir / "RefT1_vs_Case" / "details" / "samples_used.csv")
    assert set(used["group"].unique()) == {"RefT1", "Case"}

    report = (outdir / "REPORT.md").read_text(encoding="utf-8")
    assert "No analysis was executed" not in report


def test_ref_case_naming_an_absent_group_raises(tmp_path: Path) -> None:
    cfg = _write_config(
        tmp_path,
        ANALYSIS_MODE="phenotype",
        GROUP_DEFINITION={
            **THREE_GROUPS,
            "Ghost": {"group_label": "does_not_exist"},
        },
        PHENO_GROUP_REF="Ghost",
        PHENO_GROUP_CASE="Case",
    )
    with pytest.raises(ValueError, match="no samples"):
        run(LeanConfig.from_yaml(cfg))


def test_report_flags_a_run_that_produced_nothing() -> None:
    """A report with no executed analysis must say so, not claim success."""
    from phenorewire.reporting import export_final_report
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        outdir = Path(tmp)
        export_final_report(
            {
                "tool_name": "PhenoRewire",
                "analysis_mode": "auto",
                "outdir": str(outdir),
                "deliverables": [str(outdir / "REPORT.md")],
                "phenotype": {"executed": False},
                "temporal": {"executed": False},
            },
            outdir,
        )
        text = (outdir / "REPORT.md").read_text(encoding="utf-8")

    assert "No analysis was executed" in text
    assert "Run completed successfully" not in text
    # Must not point the reader at a Cytoscape file that was never written.
    assert "rewiring_network.graphml" not in text
