import json
from pathlib import Path

from phenorewire.reporting import export_final_report


def test_final_report_contains_summary_and_warnings(tmp_path: Path) -> None:
    report = {
        "tool_name": "PhenoRewire",
        "analysis_mode": "phenotype",
        "outdir": str(tmp_path),
        "n_samples_used": 10,
        "n_features_after_preprocessing": 20,
        "deliverables": ["a.csv"],
        "phenotype": {
            "executed": True,
            "ref_group": "A",
            "case_group": "B",
            "n_selected_features": 0,
        },
        "temporal": {"executed": False},
    }
    export_final_report(report, tmp_path)
    md = (tmp_path / "report" / "final_report.md").read_text(encoding="utf-8")
    payload = json.loads((tmp_path / "report" / "final_report.json").read_text(encoding="utf-8"))
    assert "Executive Summary" in md
    assert "Warnings" in md
    assert payload["tool_name"] == "PhenoRewire"
