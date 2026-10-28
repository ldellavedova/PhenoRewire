"""Documentation that users copy from must actually work."""
from __future__ import annotations

import re
from pathlib import Path

import yaml

from phenorewire.config import LeanConfig

REPO = Path(__file__).resolve().parents[3]


def _yaml_blocks(path: Path) -> list[str]:
    return re.findall(r"```yaml\n(.*?)```", path.read_text(encoding="utf-8"), flags=re.S)


def test_minimal_config_in_docs_validates(tmp_path):
    block = _yaml_blocks(REPO / "docs" / "configuration.md")[0]
    data = yaml.safe_load(block)
    data["OUTDIR"] = str(tmp_path / "out")
    for key in ("DATA_MATRIX", "METADATA"):
        data[key] = str(REPO / data[key])
    cfg = LeanConfig(**data)
    assert cfg.ANALYSIS_MODE == "phenotype"


def test_every_config_key_is_documented():
    text = (REPO / "docs" / "configuration.md").read_text(encoding="utf-8")
    missing = [k for k in LeanConfig.model_fields if f"`{k}" not in text]
    assert not missing, f"undocumented config keys: {missing}"


def test_bundled_configs_validate(tmp_path):
    for name in ("config_phenotype.yaml", "config_temporal.yaml"):
        cfg = LeanConfig.from_yaml(REPO / name)
        assert cfg.DATA_MATRIX.exists()


USER_DOCS = [REPO / "README.md", REPO / "data" / "README.md", *sorted((REPO / "docs").rglob("*.md"))]

# Output paths of 0.1.x; the user docs must describe the current layout only.
OLD_OUTPUT_PATHS = [
    "triage/", "final_report", "priority_rewired_nodes", "rewiring_node_null", "rewiring_null.csv",
    "phenotype_selection/", "rewiring_pheno", "selected_pheno_features", "node_metrics.csv",
    "X_corr_matrix", "metadata_aligned", "2dnetwork/",
]


def test_user_docs_describe_the_current_output_layout():
    stale = [
        f"{doc.relative_to(REPO)}: {old}"
        for doc in USER_DOCS
        for old in OLD_OUTPUT_PATHS
        if old in doc.read_text(encoding="utf-8")
    ]
    assert not stale, f"docs still refer to 0.1.x output paths: {stale}"


def test_relative_links_in_docs_resolve():
    broken = []
    for doc in USER_DOCS:
        for target in re.findall(r"\]\(([^)#\s]+)(?:#[^)]*)?\)", doc.read_text(encoding="utf-8")):
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            if not (doc.parent / target).exists():
                broken.append(f"{doc.relative_to(REPO)} -> {target}")
    assert not broken, f"broken links: {broken}"
