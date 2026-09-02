"""
Command-line interface for PhenoRewire.

Entry point: phenorewire --config config.yaml [options]
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phenorewire",
        description="PhenoRewire: phenotype-aware metabolic network rewiring and dynamics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  phenorewire --config config_phenotype.yaml
  phenorewire --config config_phenotype.yaml --run-2dnetwork
  phenorewire --config config_phenotype.yaml --plot-networks
  phenorewire --config config.yaml --log-level DEBUG
""",
    )
    p.add_argument("--config", required=True, help="Path to YAML configuration file.")
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    p.add_argument(
        "--run-2dnetwork",
        action="store_true",
        default=False,
        help=(
            "After the main pipeline, run the 2D network annotation expansion module. "
            "Requires SPECREBOOT_NETWORK in the config. "
            "Uses PHENOREWIRE_NETWORK if set, otherwise auto-detects the rewiring GraphML "
            "produced by this run (rewiring_network_phenotype.graphml or rewiring_network_temporal.graphml)."
        ),
    )
    p.add_argument(
        "--plot-networks",
        action="store_true",
        default=False,
        help="Produce publication-quality PNG and interactive HTML network plots.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    from .config import LeanConfig
    from .run import run, _setup_logging

    _setup_logging(args.log_level)

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        print(f"ERROR: Config file not found: {cfg_path}", file=sys.stderr)
        return 1

    config = LeanConfig.from_yaml(cfg_path)
    run(config)

    if args.run_2dnetwork or args.plot_networks:
        _run_extras(config, args, cfg_path)

    return 0


def _run_extras(config, args, cfg_path: Path) -> None:
    logger = logging.getLogger("phenorewire.cli")
    outdir = Path(config.OUTDIR)

    if args.plot_networks:
        try:
            from .twodn.viz import plot_rewiring_network
        except ImportError:
            logger.error(
                "--plot-networks requires matplotlib and pyvis. "
                "Install with: pip install 'phenorewire[viz]'"
            )
            return

        _mode = str(getattr(config, "ANALYSIS_MODE", "phenotype")).lower()
        _graphml_candidates = (
            [outdir / "triage" / "rewiring_network_temporal.graphml"]
            if _mode == "temporal"
            else [
                outdir / "triage" / "rewiring_network_phenotype.graphml",
                outdir / "triage" / "rewiring_network.graphml",  # legacy fallback
            ]
        )
        rewiring_graphml = next((p for p in _graphml_candidates if p.exists()), None)
        if rewiring_graphml is None:
            logger.warning("--plot-networks: rewiring_network graphml not found in %s/triage/", outdir)
        else:
            rewiring_graphml = Path(rewiring_graphml)
            plot_outdir = outdir / "plots"
            plot_outdir.mkdir(parents=True, exist_ok=True)
            try:
                priority_csv = outdir / "triage" / "priority_rewired_nodes.csv"
                triage_df = pd.read_csv(priority_csv) if priority_csv.exists() else pd.DataFrame()
                rewiring_df = triage_df
                plot_rewiring_network(
                    rewiring_graphml,
                    triage_df=triage_df,
                    rewiring_df=rewiring_df,
                    output_path=plot_outdir / "rewiring_network",
                )
                logger.info("Network plots written to %s", plot_outdir)
            except Exception as exc:
                logger.error("--plot-networks failed: %s", exc)

    if args.run_2dnetwork:
        try:
            from .twodn.chemical_integration import build_annotation_expansion
        except ImportError:
            logger.error(
                "--run-2dnetwork could not import the 2D annotation module. "
                "Reinstall PhenoRewire in your active environment with: pip install -e ."
            )
            return

        specreboot_path = getattr(config, "SPECREBOOT_NETWORK", None)
        phenorewire_path = getattr(config, "PHENOREWIRE_NETWORK", None)

        if specreboot_path is None:
            logger.error(
                "--run-2dnetwork: SPECREBOOT_NETWORK is not set in your config. "
                "Add 'SPECREBOOT_NETWORK: path/to/specreboot_network.graphml' to %s.",
                cfg_path,
            )
            return

        if phenorewire_path is None:
            _mode = str(getattr(config, "ANALYSIS_MODE", "phenotype")).lower()
            _candidates = (
                [outdir / "triage" / "rewiring_network_temporal.graphml"]
                if _mode == "temporal"
                else [
                    outdir / "triage" / "rewiring_network_phenotype.graphml",
                    outdir / "triage" / "rewiring_network.graphml",  # legacy fallback
                ]
            )
            _found = next((str(p) for p in _candidates if p.exists()), None)
            if _found is None:
                logger.error(
                    "--run-2dnetwork: PHENOREWIRE_NETWORK not set and no rewiring network "
                    "graphml found in %s/triage/. Run the main pipeline first.",
                    outdir,
                )
                return
            phenorewire_path = _found

        priority_csv = outdir / "triage" / "priority_rewired_nodes.csv"
        priority_df = pd.read_csv(priority_csv) if priority_csv.exists() else pd.DataFrame()

        priority_features_df = priority_df.copy()
        triage_features_csv = next(
            (p for p in outdir.rglob("priority_features.csv") if "triage" in str(p)), None
        )
        if priority_features_df.empty and triage_features_csv:
            priority_features_df = pd.read_csv(triage_features_csv)
            logger.warning(
                "--run-2dnetwork: priority_rewired_nodes.csv not found; "
                "falling back to %s for feature ranking.",
                triage_features_csv,
            )

        ann_outdir = outdir / "2dnetwork"
        ann_outdir.mkdir(parents=True, exist_ok=True)

        top_n = int(getattr(config, "ANNOTATION_TOP_N", 20))
        hop_depth = int(getattr(config, "MN_HOP_DEPTH", 2))
        min_cosine = float(getattr(config, "MN_MIN_COSINE", 0.70))

        try:
            summary_df, expansion_df, first_neighbours_df = build_annotation_expansion(
                phenorewire_graphml=Path(str(phenorewire_path)),
                specreboot_graphml=Path(str(specreboot_path)),
                priority_features_df=priority_features_df,
                rewiring_df=priority_df if not priority_df.empty else priority_features_df,
                top_n=top_n,
                hop_depth=hop_depth,
                min_cosine=min_cosine,
                output_dir=ann_outdir,
            )
            logger.info(
                "2D network annotation expansion complete: %d features, %d expansion rows, %d first-neighbour rows. "
                "Outputs in %s",
                len(summary_df),
                len(expansion_df),
                len(first_neighbours_df),
                ann_outdir,
            )
        except Exception as exc:
            logger.error("--run-2dnetwork failed: %s", exc)


if __name__ == "__main__":
    raise SystemExit(main())
