"""The two bundled examples run as documented and find their known answer.

Both datasets are synthetic (dev/benchmarks/fixtures/generate_example_data.py):
F001-F010 change in abundance and rewire, F011-F020 change in abundance only,
F021-F030 rewire without changing in abundance, F031-F060 are background.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phenorewire.config import LeanConfig
from phenorewire.run import run

REPO = Path(__file__).resolve().parents[3]

REWIRED = {f"F{i:03d}" for i in list(range(1, 11)) + list(range(21, 31))}
ABUNDANCE_ONLY = {f"F{i:03d}" for i in range(11, 21)}


@pytest.mark.parametrize(
    "config_name, comparison",
    [("config_phenotype.yaml", "Reference_vs_Case"), ("config_temporal.yaml", "Exposed_0h_vs_10h")],
)
def test_bundled_example_finds_its_known_answer(tmp_path, config_name, comparison):
    config = LeanConfig.from_yaml(REPO / config_name)
    run(config.model_copy(update={"OUTDIR": tmp_path / "out"}))

    comp = tmp_path / "out" / comparison
    assert pd.read_csv(comp / "rewiring_vs_chance.csv").loc[0, "p_empirical"] <= 0.05

    rewired = pd.read_csv(comp / "rewired_metabolites.csv")
    called = set(rewired.loc[rewired["more_than_chance"], "feature_id"])
    assert len(called & REWIRED) >= 15, sorted(called)
    assert not called - REWIRED, f"called but not rewired in the design: {sorted(called - REWIRED)}"

    routes = dict(zip(rewired["feature_id"], rewired["selection_route"]))
    assert all(routes[f] in {"connectivity", "both"} for f in called if f >= "F021")
    assert not ABUNDANCE_ONLY & called
