"""What a user sees in the terminal: readable errors, the input check, the run."""
from __future__ import annotations

from pathlib import Path

import yaml

from phenorewire.cli import main
from phenorewire.console import fmt_p
from phenorewire.layout import phenotype_comparison, safe_name, temporal_comparison

REPO = Path(__file__).resolve().parents[3]


def _config(tmp_path: Path, **overrides) -> Path:
    data = yaml.safe_load((REPO / "config_phenotype.yaml").read_text(encoding="utf-8"))
    for key in ("DATA_MATRIX", "METADATA", "SPECREBOOT_NETWORK"):
        data[key] = str(REPO / data[key])
    data["OUTDIR"] = str(tmp_path / "out")
    data.update(overrides)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


def test_missing_input_file_is_one_readable_line(tmp_path, capsys):
    cfg = _config(tmp_path, DATA_MATRIX="no_such_table.csv")
    assert main(["--config", str(cfg)]) == 1
    err = capsys.readouterr().err
    assert "Traceback" not in err
    assert "DATA_MATRIX: file not found" in err
    assert "no_such_table.csv" in err


def test_unknown_config_key_is_named(tmp_path, capsys):
    cfg = _config(tmp_path, LOG_TRANSFORMS=True)
    assert main(["--config", str(cfg)]) == 1
    err = capsys.readouterr().err
    assert "LOG_TRANSFORMS" in err and "Traceback" not in err


def test_check_reports_and_writes_nothing(tmp_path, capsys):
    cfg = _config(tmp_path)
    assert main(["--config", str(cfg), "--check"]) == 0
    out = capsys.readouterr().out
    assert "60 features x 40 samples" in out
    assert "Result: ready to run." in out
    assert not (tmp_path / "out").exists()


def test_check_fails_with_exit_code_on_bad_data(tmp_path, capsys):
    cfg = _config(tmp_path, META_SAMPLE_COL="filename")
    assert main(["--config", str(cfg), "--check"]) == 1
    assert "no column 'filename'" in capsys.readouterr().out


def test_folder_names_are_safe():
    assert safe_name("Case / treated") == "Case_treated"
    assert phenotype_comparison(Path("o"), "Ref", "Case").root == Path("o/Ref_vs_Case")
    assert temporal_comparison(Path("o"), "Cohort", "0h", "10h").root == Path("o/Cohort_0h_vs_10h")


def test_p_values_read_like_p_values():
    assert fmt_p(0.004975124378109453) == "0.005"
    assert fmt_p(0.0001) == "< 0.001"
    assert fmt_p(None) == "NA"
