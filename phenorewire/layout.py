# ── PhenoRewire · Output layout ──────────────────────────────────────────────
# The one place that decides where every output file goes.  Everything a user
# opens first sits at the top of a comparison folder; everything else is under
# details/ or run_info/.
#
#   OUTDIR/
#     REPORT.md                      start here
#     report.json                    the same, machine-readable
#     <comparison>/                  e.g. Reference_vs_Case/ or Treated_0h_vs_10h/
#       rewiring_vs_chance.csv       is the rewiring more than an arbitrary split?
#       rewired_metabolites.csv      one row per metabolite, strongest evidence first
#       selected_metabolites.csv     which metabolites entered the networks, and why
#       rewiring_network.graphml     both networks in one file, for Cytoscape
#       details/
#         selection/                 every feature's test statistics
#         network_<state>/           nodes, edges, communities, null model
#         rewiring/                  every edge with its state, global counts
#         samples_used.csv
#     run_info/                      log, input check, config, preprocessed data
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

REPORT_MD = "REPORT.md"
REPORT_JSON = "report.json"


def safe_name(label: str) -> str:
    """A group or timepoint label made safe for a folder name."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(label).strip()).strip("_")
    return s or "group"


@dataclass(frozen=True)
class RunInfo:
    root: Path

    @property
    def log(self) -> Path:
        return self.root / "run.log"

    @property
    def data_check(self) -> Path:
        return self.root / "data_check.txt"

    @property
    def config(self) -> Path:
        return self.root / "config_used.yaml"


@dataclass(frozen=True)
class Comparison:
    """Paths of one two-state comparison (phenotype groups or timepoints)."""

    root: Path

    # What the user opens
    @property
    def rewiring_vs_chance(self) -> Path:
        return self.root / "rewiring_vs_chance.csv"

    @property
    def rewired_metabolites(self) -> Path:
        return self.root / "rewired_metabolites.csv"

    @property
    def selected_metabolites(self) -> Path:
        return self.root / "selected_metabolites.csv"

    @property
    def rewiring_network(self) -> Path:
        return self.root / "rewiring_network.graphml"

    # Supporting detail
    @property
    def details(self) -> Path:
        return self.root / "details"

    @property
    def selection(self) -> Path:
        return self.details / "selection"

    @property
    def rewiring(self) -> Path:
        return self.details / "rewiring"

    @property
    def samples_used(self) -> Path:
        return self.details / "samples_used.csv"

    def network(self, state: str) -> Path:
        return self.details / f"network_{safe_name(state)}"

    # Optional extras (--run-2dnetwork, --plot-networks)
    @property
    def molecular_network(self) -> Path:
        return self.root / "molecular_network"

    @property
    def plots(self) -> Path:
        return self.root / "plots"


def run_info(outdir: Path) -> RunInfo:
    return RunInfo(Path(outdir) / "run_info")


def phenotype_comparison(outdir: Path, ref_group: str, case_group: str) -> Comparison:
    return Comparison(Path(outdir) / f"{safe_name(ref_group)}_vs_{safe_name(case_group)}")


def temporal_comparison(outdir: Path, group: str, t1_label: str, t2_label: str) -> Comparison:
    return Comparison(Path(outdir) / f"{safe_name(group)}_{safe_name(t1_label)}_vs_{safe_name(t2_label)}")


def find_comparisons(outdir: Path) -> list[Comparison]:
    """Comparison folders of a finished run, found by their rewiring network."""
    return [
        Comparison(p.parent)
        for p in sorted(Path(outdir).glob("*/rewiring_network.graphml"))
    ]
