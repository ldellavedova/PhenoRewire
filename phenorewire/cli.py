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
  phenorewire --config config_phenotype.yaml --check     check the input data only
  phenorewire --config config_phenotype.yaml             run the analysis
  phenorewire --config config_phenotype.yaml --run-2dnetwork
  phenorewire --config config_phenotype.yaml --plot-networks

Results: OUTDIR/REPORT.md says what was found and which file to open first.
""",
    )
    p.add_argument("--config", required=True, help="Path to YAML configuration file.")
    p.add_argument(
        "--check",
        action="store_true",
        default=False,
        help="Only check the input data (samples, groups, timepoints, values) and stop.",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Console verbosity (default: INFO). The full log is always in OUTDIR/run_info/run.log.",
    )
    p.add_argument(
        "--run-2dnetwork",
        action="store_true",
        default=False,
        help=(
            "After the main pipeline, run the 2D network annotation expansion module. "
            "Requires SPECREBOOT_NETWORK in the config. "
            "Uses PHENOREWIRE_NETWORK if set, otherwise the rewiring network of each comparison of this run."
        ),
    )
    p.add_argument(
        "--plot-networks",
        action="store_true",
        default=False,
        help="Produce publication-quality PNG and interactive HTML network plots.",
    )
    return p


def _config_error_text(exc: Exception) -> str:
    """A pydantic validation error as one line per problem."""
    errors = getattr(exc, "errors", None)
    if not callable(errors):
        return str(exc)
    lines = []
    for err in errors():
        where = ".".join(str(x) for x in err.get("loc", ())) or "config"
        msg = str(err.get("msg", "")).removeprefix("Value error, ")
        ctx = err.get("ctx") or {}
        if "error" in ctx and str(ctx["error"]) not in msg:
            msg = f"{msg}: {ctx['error']}"
        lines.append(f"- {where}: {msg}")
    return "\n".join(lines)


def _fail(title: str, detail: str, debug: bool) -> int:
    print(f"\nERROR: {title}", file=sys.stderr)
    if detail:
        print(detail, file=sys.stderr)
    if not debug:
        print("\n(Run with --log-level DEBUG to see the full traceback.)", file=sys.stderr)
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    from pydantic import ValidationError

    from .checks import check_inputs
    from .config import LeanConfig
    from .console import setup_console
    from .run import run

    debug = args.log_level == "DEBUG"
    setup_console(args.log_level)

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        return _fail(f"config file not found: {cfg_path}", "", debug=True)

    try:
        config = LeanConfig.from_yaml(cfg_path)
    except ValidationError as exc:
        if debug:
            raise
        return _fail(
            f"the config file {cfg_path} has a problem. Relative paths are read from the "
            "config file's own folder.",
            _config_error_text(exc),
            debug,
        )
    except (ValueError, OSError) as exc:
        if debug:
            raise
        return _fail(f"could not read the config file {cfg_path}.", str(exc), debug)

    if args.check:
        chk = check_inputs(config)
        print(chk.text())
        return 0 if chk.ok else 1

    try:
        run(config, config_path=cfg_path)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        if debug:
            raise
        return _fail(
            "the run stopped.",
            f"{exc}\n\nThe log of the run so far: {Path(config.OUTDIR) / 'run_info' / 'run.log'}",
            debug,
        )

    if args.run_2dnetwork or args.plot_networks:
        _run_extras(config, args, cfg_path)

    return 0


def _run_extras(config, args, cfg_path: Path) -> None:
    from .console import rel
    from .layout import find_comparisons

    logger = logging.getLogger("phenorewire.cli")
    outdir = Path(config.OUTDIR)
    comparisons = find_comparisons(outdir)
    if not comparisons:
        logger.warning("No rewiring network found in %s; nothing to annotate or plot.", outdir)
        return

    if args.plot_networks:
        try:
            from .twodn.viz import plot_rewiring_network
        except ImportError:
            logger.error(
                "--plot-networks requires matplotlib and pyvis. "
                "Install with: pip install 'phenorewire[viz]'"
            )
            return

        for comp in comparisons:
            comp.plots.mkdir(parents=True, exist_ok=True)
            try:
                ranked = pd.read_csv(comp.rewired_metabolites) if comp.rewired_metabolites.exists() else pd.DataFrame()
                plot_rewiring_network(
                    comp.rewiring_network,
                    triage_df=ranked,
                    rewiring_df=ranked,
                    output_path=comp.plots / "rewiring_network",
                )
                logging.getLogger("phenorewire.progress").info("  Network plots: %s/", rel(comp.plots))
            except Exception as exc:
                logger.error("--plot-networks failed for %s: %s", comp.root.name, exc)

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

        top_n = int(getattr(config, "ANNOTATION_TOP_N", 20))
        hop_depth = int(getattr(config, "MN_HOP_DEPTH", 2))
        min_cosine = float(getattr(config, "MN_MIN_COSINE", 0.70))

        for comp in comparisons:
            ranked = pd.read_csv(comp.rewired_metabolites) if comp.rewired_metabolites.exists() else pd.DataFrame()
            graphml = Path(str(phenorewire_path)) if phenorewire_path else comp.rewiring_network
            ann_outdir = comp.molecular_network
            ann_outdir.mkdir(parents=True, exist_ok=True)
            try:
                summary_df, expansion_df, first_neighbours_df = build_annotation_expansion(
                    phenorewire_graphml=graphml,
                    specreboot_graphml=Path(str(specreboot_path)),
                    priority_features_df=ranked,
                    rewiring_df=ranked,
                    top_n=top_n,
                    hop_depth=hop_depth,
                    min_cosine=min_cosine,
                    output_dir=ann_outdir,
                )
                logging.getLogger("phenorewire.progress").info(
                    "  Molecular-network annotation: %d metabolites, %d neighbours -> %s/",
                    len(summary_df), len(first_neighbours_df), rel(ann_outdir),
                )
            except Exception as exc:
                logger.error("--run-2dnetwork failed for %s: %s", comp.root.name, exc)


if __name__ == "__main__":
    raise SystemExit(main())
