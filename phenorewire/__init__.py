# PhenoRewire — interpretable metabolomics rewiring pipeline
#
# Pipeline execution order:
#   config.py              →  validate YAML config (Pydantic v2)
#   preprocessing.py       →  Step 1: filter, normalize (NORM_METHOD), log2 transform
#   phenotype_selection.py →  Step 2a: feature selection by abundance (phenotype mode)
#   temporal_selection.py  →  Step 2b: feature selection by abundance (temporal mode)
#   connectivity.py        →  Step 2c: selection by change in co-variation partners
#   networks.py            →  Step 3: network construction + rewiring quantification
#   nulls.py               →  Step 4: rewiring null, edge null, modularity vs configuration
#   triage.py              →  Step 5: per-metabolite tables, ranked
#   reporting.py           →  Step 6: REPORT.md + report.json
#
# Around them: checks.py (input check), console.py (what the terminal shows),
# layout.py (where every output file goes), stats.py and utils.py (shared
# statistics and helpers).
#
# Entry point: run.py (pipeline orchestrator) — invoked via CLI (cli.py)
# CLI command: phenorewire --config your_config.yaml

from ._version import __version__
from .run import run
from .cli import main
from .config import LeanConfig

__all__ = [
    "run",
    "main",
    "LeanConfig",
    "__version__",
]
