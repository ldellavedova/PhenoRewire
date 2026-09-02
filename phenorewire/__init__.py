# PhenoRewire — interpretable metabolomics rewiring pipeline
#
# Pipeline execution order:
#   config.py              →  validate YAML config (Pydantic v2)
#   preprocessing.py       →  Step 1: filter, median-normalize, log2 transform
#   phenotype_selection.py →  Step 2a: feature selection (phenotype mode)
#   temporal_selection.py  →  Step 2b: feature selection (temporal mode)
#   networks.py            →  Step 3: network construction + rewiring quantification
#   triage.py              →  Step 4: triage scoring and priority output
#   reporting.py           →  Step 5: final report (Markdown + JSON)
#
# Entry point: run.py (pipeline orchestrator) — invoked via CLI (cli.py)
# CLI command: phenorewire --config your_config.yaml

from .run import run
from .cli import main
from .config import LeanConfig

__version__ = "0.1.1"

__all__ = [
    "run",
    "main",
    "LeanConfig",
    "__version__",
]
