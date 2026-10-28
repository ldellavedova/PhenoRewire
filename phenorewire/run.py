from __future__ import annotations

from pathlib import Path
from typing import Optional
import logging
import time

import numpy as np
import pandas as pd
import yaml

from . import layout
from ._version import __version__
from .checks import (  # noqa: F401  (operating range re-exported for callers of run.py)
    OPERATING_RANGE_MIN_N,
    OPERATING_RANGE_RECOMMENDED_N,
    InputDataError,
    check_inputs,
    operating_range,
)
from .config import LeanConfig
from .console import Steps, add_file_log, has_any_handler, progress, rel, remove_handler, setup_console
from .preprocessing import prepare_matrices
from .phenotype_selection import phenotype_selection
from .networks import (
    correlation_network,
    export_network,
    compute_network_metrics,
    compute_rewiring,
    summarize_rewiring,
    export_rewiring_network,
    check_louvain_stability,
)
from .connectivity import connectivity_selection
from .nulls import modularity_null, network_edge_null, rewiring_label_null
from .temporal_selection import temporal_selection
from .triage import rank_network_nodes, rewired_metabolites
from .reporting import comparison_verdict, export_final_report
from .utils import NORM_LABELS, ensure_dir, parse_time_numeric, timepoint_label

logger = logging.getLogger(__name__)

# Selection permutations are drawn from one generator consumed feature by
# feature (stats.spearman_permutation_test, stream="shared").  This is the
# scheme every published PhenoRewire analysis was run with, so abundance-based
# selection reproduces those analyses exactly: v0.2.0 keeps every feature they
# selected and only adds features through connectivity-based selection.
PERMUTATION_STREAM = "shared"


# ============================================================
# Logging
# ============================================================

def _setup_logging(level: str = "INFO") -> None:
    setup_console(level)


# ============================================================
# Helpers
# ============================================================

def _subset_samples(meta2: pd.DataFrame, sample_col: str, mask: pd.Series) -> list[str]:
    return meta2.loc[mask, sample_col].astype(str).tolist()


def _safe_mean_rows(X: pd.DataFrame) -> pd.Series:
    if X is None or X.empty:
        return pd.Series(dtype=float)
    return X.mean(axis=1, skipna=True)


def _add_directional_columns(
    selected: pd.DataFrame,
    X_corr: pd.DataFrame,
    meta2: pd.DataFrame,
    *,
    sample_col: str,
    ref_group: str,
    case_group: str,
) -> pd.DataFrame:
    sel = selected.copy()
    feats = sel.index.astype(str).tolist()
    Xs = X_corr.loc[feats, :]

    samp_ref = _subset_samples(meta2, sample_col, meta2["group"] == ref_group)
    samp_case = _subset_samples(meta2, sample_col, meta2["group"] == case_group)

    mean_ref = _safe_mean_rows(Xs.loc[:, samp_ref]) if len(samp_ref) else pd.Series(0.0, index=feats)
    mean_case = _safe_mean_rows(Xs.loc[:, samp_case]) if len(samp_case) else pd.Series(0.0, index=feats)

    sel["mean_ref"] = mean_ref.reindex(feats).values
    sel["mean_case"] = mean_case.reindex(feats).values
    sel["log2FC_case_vs_ref"] = sel["mean_case"].astype(float) - sel["mean_ref"].astype(float)

    eps = 1e-12
    sel["direction"] = np.where(
        sel["log2FC_case_vs_ref"].astype(float) > eps,
        f"up_in_{case_group}",
        np.where(sel["log2FC_case_vs_ref"].astype(float) < -eps, f"up_in_{ref_group}", "flat"),
    )
    return sel


def _check_group_sample_sizes(
    group_series: "pd.Series",
    *,
    min_hard_stop: int = 5,
    min_warn: int = 10,
    max_imbalance_ratio: float = 3.0,
) -> None:
    """Validate per-group sample counts; warn or error as appropriate."""
    counts = group_series.value_counts()
    for grp, n in counts.items():
        if n < min_hard_stop:
            raise ValueError(
                f"Group '{grp}' has only {n} sample(s). At least {min_hard_stop} samples per "
                "group are required for reliable Spearman permutation testing. "
                "Options: add replicates, merge timepoints, or use a looser GROUP_DEFINITION."
            )
        if n < min_warn:
            logger.warning(
                "Group '%s' has %d samples (fewer than the recommended %d). "
                "Permutation p-values and network topology will be unstable; "
                "interpret results with caution.",
                grp, n, min_warn,
            )

    if len(counts) >= 2:
        max_n = int(counts.max())
        min_n = int(counts.min())
        if min_n > 0 and max_n / min_n > max_imbalance_ratio:
            logger.warning(
                "Group size imbalance: largest group (%d samples) is %.1fx the smallest (%d). "
                "Spearman correlations with the binary group label may be biased toward the "
                "larger group.  Consider balanced subsampling if imbalance is >3:1.",
                max_n, max_n / min_n, min_n,
            )


def _log_operating_range(label: str, op: dict) -> None:
    # The input check has already told the user; the log keeps the record.
    if op["status"] == "below":
        logger.debug(
            "%s: %d samples in the smaller state, below the operating range (%d).",
            label, op["min_n_per_state"], OPERATING_RANGE_MIN_N,
        )
    elif op["status"] == "reduced":
        logger.debug(
            "%s: %d samples in the smaller state (recommended >= %d).",
            label, op["min_n_per_state"], OPERATING_RANGE_RECOMMENDED_N,
        )


def _apply_connectivity_selection(
    config: LeanConfig,
    abundance_all: pd.DataFrame,
    abundance_selected: pd.DataFrame,
    X_a: pd.DataFrame,
    X_b: pd.DataFrame,
    outdir: Path,
    *,
    q_col: str,
    label: str,
) -> tuple[pd.DataFrame, dict]:
    """Add connectivity-based selection to the abundance-based one, when the data allow.

    ``X_a`` / ``X_b`` hold every preprocessed feature for the two states that
    will be compared.  Connectivity-based selection runs automatically when both
    states have at least ``CONNECTIVITY.min_n_per_state`` samples (and the run is
    not reproducing an earlier version).  Otherwise the abundance selection is
    returned untouched and the reason is logged and reported.

    When it runs, the returned table holds the features selected by either
    route, each with ``selection_route`` (abundance / connectivity / both), the
    connectivity statistics, and ``q_selection``, the smaller of the two
    q-values, which the triage score uses.
    """
    cc = config.CONNECTIVITY
    n_min = int(min(X_a.shape[1], X_b.shape[1]))
    meta: dict = {"min_n_per_state": n_min, "required_n": int(cc.min_n_per_state)}

    if n_min < cc.min_n_per_state:
        msg = (
            f"{label}: features were selected by abundance only. Selecting by change in "
            f"co-variation partners needs at least {cc.min_n_per_state} samples in each state, "
            f"and the smaller state has {n_min}; below that it finds almost nothing in "
            "simulation and its permutation null is too coarse to trust."
        )
        logger.warning("%s", msg)
        meta.update({"connectivity_run": False, "reason": msg})
        if len(abundance_selected) < int(config.MIN_FEATURES_HARD_STOP):
            raise ValueError(
                f"{label}: only {len(abundance_selected)} feature(s) selected, fewer than "
                f"MIN_FEATURES_HARD_STOP = {config.MIN_FEATURES_HARD_STOP}."
            )
        return abundance_selected, meta

    abundance_ids = set(abundance_selected.index.astype(str))
    conn_table, conn_summary = connectivity_selection(
        X_a, X_b,
        top_k=cc.top_k,
        n_permutations=cc.n_permutations,
        fdr_alpha=cc.fdr_alpha,
        seed=config.RANDOM_SEED,
    )
    outdir.mkdir(parents=True, exist_ok=True)
    conn_table.join(abundance_all[[c for c in ("name", q_col) if c in abundance_all.columns]], how="left") \
        .to_csv(outdir / "all_features_connectivity.csv")
    connectivity_ids = set(conn_table.index[conn_table["selected_connectivity"]])
    meta.update({"connectivity_run": True, **conn_summary})

    chosen = abundance_ids | connectivity_ids
    order = [f for f in abundance_all.index.astype(str) if f in chosen]
    out = abundance_all.loc[order].copy()
    out["selection_route"] = [
        "both" if (f in abundance_ids and f in connectivity_ids)
        else "abundance" if f in abundance_ids else "connectivity"
        for f in order
    ]
    q_abund = out[q_col].astype(float) if q_col in out.columns else pd.Series(1.0, index=out.index)
    out = out.join(conn_table, how="left")
    q_conn = out["q_connectivity"].astype(float)
    out["q_selection"] = np.where(
        out["selection_route"] == "abundance", q_abund,
        np.where(out["selection_route"] == "connectivity", q_conn, np.fmin(q_abund, q_conn)),
    )

    counts = out["selection_route"].value_counts()
    meta.update({
        "n_selected": int(len(out)),
        "n_abundance_only": int(counts.get("abundance", 0)),
        "n_connectivity_only": int(counts.get("connectivity", 0)),
        "n_both": int(counts.get("both", 0)),
    })
    logger.debug(
        "%s: %d features selected (%d abundance only, %d connectivity only, %d both).",
        label, meta["n_selected"], meta["n_abundance_only"], meta["n_connectivity_only"], meta["n_both"],
    )
    if len(out) < int(config.MIN_FEATURES_HARD_STOP):
        raise ValueError(
            f"{label}: only {len(out)} feature(s) selected by abundance and connectivity, fewer than "
            f"MIN_FEATURES_HARD_STOP = {config.MIN_FEATURES_HARD_STOP}."
        )
    return out, meta


def _selection_line(n_selected: int, meta: dict, what: str) -> str:
    if meta.get("connectivity_run"):
        return (
            f"{n_selected} metabolites selected: {meta.get('n_abundance_only', 0)} by {what} only, "
            f"{meta.get('n_connectivity_only', 0)} by change in co-variation partners only, "
            f"{meta.get('n_both', 0)} by both"
        )
    return f"{n_selected} metabolites selected by {what}"


def _build_and_export_network(
    selected_features: pd.DataFrame,
    X_subset: pd.DataFrame,
    outdir: Path,
    *,
    abs_r_thr: float,
    fdr_alpha: float,
    adaptive_thresholds: bool,
    min_edges: int,
    max_density: float,
    floor_abs_r_thr: float,
    ceil_fdr_alpha: float,
    top_n: int,
    seed: int = 42,
    topology_weight: float = 0.7,
    selection_weight: float = 0.3,
    louvain_stability_check: bool = False,
    null_permutations: int = 0,
    null_configuration_models: int = 0,
) -> tuple[pd.DataFrame, dict]:
    _, edges_df, threshold_meta = correlation_network(
        X_feat_by_samp=X_subset,
        abs_r_thr=abs_r_thr,
        fdr_alpha=fdr_alpha,
        adaptive_thresholds=adaptive_thresholds,
        min_edges=min_edges,
        max_density=max_density,
        floor_abs_r_thr=floor_abs_r_thr,
        ceil_fdr_alpha=ceil_fdr_alpha,
    )
    nodes_annotated, community_summary, network_summary = compute_network_metrics(selected_features, edges_df, seed=seed)
    if not network_summary.empty:
        for key, value in threshold_meta.items():
            network_summary[key] = value
    export_network(nodes_annotated, edges_df, outdir)
    community_summary.to_csv(outdir / "communities.csv", index=False)
    network_summary.to_csv(outdir / "network_summary.csv", index=False)
    # nodes.csv: every node with its metrics and triage score, best first.
    ranked = rank_network_nodes(
        nodes_annotated, topology_weight=topology_weight, selection_weight=selection_weight
    )
    if not ranked.empty:
        ranked.to_csv(outdir / "nodes.csv", index=False)

    louvain_stability_meta: dict = {}
    if louvain_stability_check:
        try:
            louvain_stability_meta = check_louvain_stability(edges_df, n_runs=10, base_seed=seed)
            if louvain_stability_meta.get("warning"):
                logger.warning("%s", louvain_stability_meta["warning"])
        except ImportError:
            logger.warning(
                "LOUVAIN_STABILITY_CHECK is enabled but scikit-learn is not installed. "
                "Stability check skipped. Install with: pip install scikit-learn"
            )
            louvain_stability_meta = {}

    summary_row = network_summary.iloc[0].to_dict() if not network_summary.empty else {}
    report = {
        "n_nodes": int(summary_row.get("n_nodes", 0)),
        "n_edges": int(summary_row.get("n_edges", 0)),
        "density": float(summary_row.get("density", 0.0)),
        "modularity": float(summary_row.get("modularity", 0.0)),
        "n_communities": int(summary_row.get("n_communities", 0)),
        "adaptive_thresholds_used": bool(summary_row.get("adaptive_thresholds_used", False)),
        "applied_abs_r_thr": float(summary_row.get("applied_abs_r_thr", abs_r_thr)),
        "applied_fdr_alpha": float(summary_row.get("applied_fdr_alpha", fdr_alpha)),
        "base_abs_r_thr": float(summary_row.get("base_abs_r_thr", abs_r_thr)),
        "base_fdr_alpha": float(summary_row.get("base_fdr_alpha", fdr_alpha)),
        "target_min_edges": int(summary_row.get("target_min_edges", 0)),
        "target_max_edges": int(summary_row.get("target_max_edges", 0)),
        "priority_features": ranked.head(top_n).to_dict(orient="records"),
        "priority_columns": ranked.columns.tolist(),
        "louvain_stability": louvain_stability_meta,
        "n_samples": int(X_subset.shape[1]),
    }

    if null_permutations > 0:
        null_meta = network_edge_null(
            X_subset,
            abs_r_thr=report["applied_abs_r_thr"],
            fdr_alpha=report["applied_fdr_alpha"],
            n_permutations=null_permutations,
            seed=seed,
        )
        if null_meta["observed_edges"] != report["n_edges"]:
            raise RuntimeError(
                f"Edge null rebuilt {null_meta['observed_edges']} edges, the network has {report['n_edges']}."
            )
        null_meta.update(modularity_null(edges_df, n_reps=null_configuration_models, seed=seed))
        pd.DataFrame([null_meta]).to_csv(outdir / "null_model.csv", index=False)
        report["null_model"] = null_meta
    return edges_df, report


def _check_network_hard_stops(
    network_report: dict,
    label: str,
    min_edges: int,
    min_nodes: int,
) -> None:
    """Raise an informative error if a network is too sparse to be meaningful."""
    n_edges = int(network_report.get("n_edges", 0))
    n_nodes = int(network_report.get("n_nodes", 0))
    applied_r = network_report.get("applied_abs_r_thr", "?")
    applied_q = network_report.get("applied_fdr_alpha", "?")

    if n_edges < min_edges:
        raise ValueError(
            f"Network '{label}' has only {n_edges} edge(s) after all threshold adjustments "
            f"(applied |r| >= {applied_r}, FDR q <= {applied_q}). "
            f"Minimum required: NETWORK_MIN_EDGES_HARD_STOP = {min_edges}.\n"
            "To resolve: add more samples per group, reduce feature count (stricter presence "
            "filter), lower CORR_ABS_THRESHOLD, raise CORR_FDR_ALPHA, or set "
            "NETWORK_MIN_EDGES_HARD_STOP to a smaller value in your config."
        )
    if n_nodes < min_nodes:
        raise ValueError(
            f"Network '{label}' has only {n_nodes} node(s). "
            f"Minimum required: NETWORK_MIN_NODES_HARD_STOP = {min_nodes}.\n"
            "To resolve: select more features (lower FDR_ALPHA or enable "
            "ADAPTIVE_SELECTION_THRESHOLDS) or reduce MIN_FEATURES_HARD_STOP."
        )


# Columns of selected_metabolites.csv, in order; the full statistics of every
# feature stay in details/selection/.
_SELECTED_COLUMNS = {
    "phenotype": ["direction", "mean_ref", "mean_case", "log2FC_case_vs_ref", "r_pheno", "q_pheno"],
    "temporal": ["time_direction", "mean_t1", "mean_t2", "delta_t2_minus_t1", "r_time", "q_time"],
}


def _write_selected_table(selected: pd.DataFrame, kind: str, path: Path) -> None:
    cols = ["name", "mz", "rt", "selection_route"] + _SELECTED_COLUMNS[kind] + [
        "connectivity_score", "q_connectivity", "q_selection",
    ]
    table = selected[[c for c in cols if c in selected.columns]].copy()
    table.index = table.index.astype(str)
    table.index.name = "feature_id"
    table.to_csv(path)


def _compare_states(
    config: LeanConfig,
    comp: layout.Comparison,
    selected: pd.DataFrame,
    X_a: pd.DataFrame,
    X_b: pd.DataFrame,
    samples: pd.DataFrame,
    *,
    kind: str,
    label_a: str,
    label_b: str,
    stop_labels: tuple[str, str],
    direction_col: str,
    steps: Steps,
) -> tuple[dict, dict, dict]:
    """Networks of the two states, their comparison, and the comparison against chance.

    Shared by the phenotype and the temporal analysis; returns the two network
    reports and the rewiring report.
    """
    top_n = int(config.REPORT_TOP_N)
    null_perm = int(config.NULL_MODEL.n_permutations) if config.NULL_MODEL.enabled else 0
    null_cm = int(config.NULL_MODEL.n_configuration_models) if config.NULL_MODEL.enabled else 0
    net_kwargs = dict(
        abs_r_thr=config.CORR_ABS_THRESHOLD,
        fdr_alpha=config.CORR_FDR_ALPHA,
        adaptive_thresholds=config.ADAPTIVE_NETWORK_THRESHOLDS,
        min_edges=config.NETWORK_MIN_EDGES,
        max_density=config.NETWORK_MAX_DENSITY,
        floor_abs_r_thr=config.NETWORK_ABS_THRESHOLD_FLOOR,
        ceil_fdr_alpha=config.NETWORK_FDR_ALPHA_CEILING,
        top_n=top_n,
        seed=config.RANDOM_SEED,
        topology_weight=config.TRIAGE_TOPOLOGY_WEIGHT,
        selection_weight=config.TRIAGE_SELECTION_WEIGHT,
        louvain_stability_check=config.LOUVAIN_STABILITY_CHECK,
        null_permutations=null_perm,
        null_configuration_models=null_cm,
    )

    steps.start(f"Building the two networks: {label_a} and {label_b}")
    edges_a, report_a = _build_and_export_network(selected, X_a, comp.network(label_a), **net_kwargs)
    _check_network_hard_stops(report_a, stop_labels[0], config.NETWORK_MIN_EDGES_HARD_STOP, config.NETWORK_MIN_NODES_HARD_STOP)
    edges_b, report_b = _build_and_export_network(selected, X_b, comp.network(label_b), **net_kwargs)
    _check_network_hard_stops(report_b, stop_labels[1], config.NETWORK_MIN_EDGES_HARD_STOP, config.NETWORK_MIN_NODES_HARD_STOP)
    for lbl, rep in ((label_a, report_a), (label_b, report_b)):
        steps.detail(
            f"{lbl}: {rep['n_nodes']} metabolites, {rep['n_edges']} correlations "
            f"(|r| >= {rep['applied_abs_r_thr']:.2f}, q <= {rep['applied_fdr_alpha']:.2f})"
        )

    steps.start(f"Comparing the networks: {label_a} vs {label_b}")
    shared, only_a, only_b, union = compute_rewiring(
        edges_a, edges_b, label_a=label_a, label_b=label_b, sign_switch_min_r=config.SIGN_SWITCH_MIN_R,
    )
    rewiring_summary, rewiring_nodes = summarize_rewiring(union, label_a=label_a, label_b=label_b)
    ensure_dir(comp.rewiring)
    union.drop(columns=["edge_id"], errors="ignore").to_csv(comp.rewiring / "edges.csv", index=False)
    rewiring_summary.to_csv(comp.rewiring / "rewiring_summary.csv", index=False)
    samples.to_csv(comp.samples_used, index=False)
    steps.detail(
        f"{len(shared)} correlations shared, {len(only_a)} only in {label_a}, {len(only_b)} only in {label_b}"
    )

    rewiring: dict = {"summary": rewiring_summary.iloc[0].to_dict() if not rewiring_summary.empty else {}}
    node_null = None
    if null_perm:
        steps.detail(f"testing against chance ({null_perm} permutations of the state labels) ...")
        summary, node_null = rewiring_label_null(
            X_a,
            X_b,
            thresholds_a=(report_a["applied_abs_r_thr"], report_a["applied_fdr_alpha"]),
            thresholds_b=(report_b["applied_abs_r_thr"], report_b["applied_fdr_alpha"]),
            n_permutations=null_perm,
            seed=config.RANDOM_SEED,
        )
        if summary["observed_state_specific_edges"] != len(only_a) + len(only_b):
            # Both counts come from the same edge filter; a mismatch is a bug.
            raise RuntimeError(
                "Rewiring null rebuilt a different observed network "
                f"({summary['observed_state_specific_edges']} vs {len(only_a) + len(only_b)} state-specific edges)."
            )
        pd.DataFrame([summary]).to_csv(comp.rewiring_vs_chance, index=False)
        summary["n_nodes_q05"] = int((node_null["q_bh"] <= 0.05).sum())
        rewiring["null"] = summary

    table = rewired_metabolites(
        selected, rewiring_nodes, node_null, label_a=label_a, label_b=label_b, direction_col=direction_col,
    )
    table.to_csv(comp.rewired_metabolites, index=False)
    _write_selected_table(selected, kind, comp.selected_metabolites)
    # Node attributes: the selection columns plus every rewiring column, q_rewiring included.
    node_rewiring = table.drop(columns=["rank"] + [c for c in table.columns if c in selected.columns])
    export_rewiring_network(selected, union, node_rewiring, comp.rewiring_network)
    rewiring["priority_nodes"] = table.head(top_n).to_dict(orient="records")
    steps.detail(comparison_verdict({"rewiring": rewiring}, label_a, label_b))
    return report_a, report_b, rewiring


def _write_config_used(config: LeanConfig, path: Path) -> None:
    data = config.model_dump(mode="json")
    with path.open("w", encoding="utf-8") as f:
        f.write("# The configuration this run used, every default filled in.\n")
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


def _elapsed(seconds: float) -> str:
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds} s"
    return f"{seconds // 60} min {seconds % 60:02d} s"


# ============================================================
# Main
# ============================================================

def run(config: LeanConfig, config_path: Optional[Path] = None) -> None:
    """Run PhenoRewire with ``config``.

    Writes everything under ``config.OUTDIR``; ``REPORT.md`` there says what to
    open first.  ``config_path`` is only shown in the console and the report.
    """
    outdir = ensure_dir(Path(config.OUTDIR))
    info = layout.run_info(outdir)
    ensure_dir(info.root)
    own_console = None if has_any_handler() else setup_console("INFO")
    file_log = add_file_log(info.log)
    try:
        _run(config, outdir, info, config_path)
    except Exception:
        logger.debug("Run stopped with an error.", exc_info=True)
        raise
    finally:
        remove_handler(file_log)
        if own_console is not None:
            remove_handler(own_console)


def _run(config: LeanConfig, outdir: Path, info: layout.RunInfo, config_path: Optional[Path]) -> None:
    started = time.time()
    progress.info("PhenoRewire %s", __version__)
    if config_path is not None:
        progress.info("  config : %s", rel(Path(config_path)))
    progress.info("  output : %s", rel(outdir))
    _write_config_used(config, info.config)

    df = pd.read_csv(config.DATA_MATRIX, low_memory=False)
    meta = pd.read_csv(config.METADATA, low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]
    meta.columns = [str(c).strip() for c in meta.columns]

    # -----------------------------
    # 1) Input check
    # -----------------------------
    check = check_inputs(config, df=df, meta=meta)
    info.data_check.write_text(check.text() + "\n", encoding="utf-8")
    # input check, preprocessing, three steps per comparison, report
    steps = Steps(total=3 + 3 * max(1, len(check.comparisons)))
    steps.start("Checking the input data")
    for line in check.lines:
        steps.detail(line)
    for warning in check.warnings:
        logger.warning("%s", warning)
    if not check.ok:
        raise InputDataError(check.error_message())

    requested_mode = str(config.ANALYSIS_MODE).lower()

    # -----------------------------
    # 2) Preprocessing -> matrices
    # -----------------------------
    steps.start("Preprocessing")
    X_corr, meta2, feat_anno = prepare_matrices(
        df=df,
        meta=meta,
        feature_id_col=config.FEATURE_ID_COL,
        mz_col=config.MZ_COL,
        rt_col=config.RT_COL,
        name_col=config.NAME_COL,
        intensity_regex=config.INTENSITY_REGEX,
        meta_sample_col=config.META_SAMPLE_COL,
        outdir=info.root,
        min_nonzero_intensity=config.MIN_NONZERO_INTENSITY,
        min_presence_global=config.MIN_PRESENCE_GLOBAL,
        norm_method=config.NORM_METHOD,
        log_transform=config.LOG_TRANSFORM,
        group_definition=config.GROUP_DEFINITION,
    )
    steps.detail(
        f"{X_corr.shape[0]} of {len(df)} features kept; {NORM_LABELS[config.NORM_METHOD]}"
        + (", log2(x + 1)" if config.LOG_TRANSFORM else ", no log transform")
    )
    logger.debug("Groups used:\n%s", meta2["group"].value_counts(dropna=False).to_string())

    # keep only allowed groups
    allowed_groups = list(config.GROUP_DEFINITION.keys())
    meta2 = meta2[meta2["group"].isin(allowed_groups)].copy()
    if meta2.empty:
        raise ValueError(
            "After filtering to GROUP_DEFINITION keys, no samples remain. "
            "Check GROUP_DEFINITION and metadata columns/values."
        )

    # align X columns to meta order
    samp_ids = meta2[config.META_SAMPLE_COL].astype(str).tolist()
    X_corr = X_corr.loc[:, samp_ids]

    unique_groups = sorted(meta2["group"].unique())

    # phenotype direction
    ref_group = getattr(config, "PHENO_GROUP_REF", None)
    case_group = getattr(config, "PHENO_GROUP_CASE", None)

    wants_pheno = requested_mode in {"auto", "phenotype", "both"}
    wants_temporal = (
        requested_mode in {"auto", "temporal", "both"}
        and bool(config.TEMPORAL_CORRELATION.enabled)
    )

    # When the config names the two groups to contrast, honour that choice instead
    # of refusing to run because a third group also has samples.
    if wants_pheno and len(unique_groups) > 2 and ref_group is not None and case_group is not None:
        contrast = [str(ref_group), str(case_group)]
        missing = [g for g in contrast if g not in unique_groups]
        if missing:
            raise ValueError(
                f"PHENO_GROUP_REF/CASE name group(s) with no samples: {missing}. "
                f"Groups with samples after filtering: {unique_groups}."
            )
        temporal_group = str(config.TEMPORAL_CORRELATION.phenotype or "")
        if wants_temporal and temporal_group and temporal_group not in contrast:
            raise ValueError(
                f"Cannot restrict the phenotype contrast to {contrast}: "
                f"TEMPORAL_CORRELATION.phenotype='{temporal_group}' lies outside it, so the "
                "temporal analysis would lose its samples.\n"
                "To resolve: run the two analyses as separate jobs, or narrow GROUP_DEFINITION "
                "so only the groups you need match samples."
            )
        dropped = [g for g in unique_groups if g not in contrast]
        logger.debug(
            "PHENO_GROUP_REF/CASE select the contrast '%s' vs '%s'; ignoring %d other group(s): %s",
            contrast[0], contrast[1], len(dropped), dropped,
        )
        meta2 = meta2[meta2["group"].isin(contrast)].copy()
        samp_ids = meta2[config.META_SAMPLE_COL].astype(str).tolist()
        X_corr = X_corr.loc[:, samp_ids]
        unique_groups = sorted(meta2["group"].unique())

    can_run_pheno = wants_pheno and len(unique_groups) == 2
    can_run_temporal = wants_temporal

    # An analysis the user asked for explicitly must never be skipped in silence:
    # a run that produces no results has to fail loudly, not exit 0 with an empty report.
    if requested_mode in {"phenotype", "both"} and not can_run_pheno:
        raise ValueError(
            f"ANALYSIS_MODE='{requested_mode}' requires exactly two groups with samples, but "
            f"{len(unique_groups)} group(s) remain after filtering: {unique_groups}.\n"
            "To resolve: set PHENO_GROUP_REF and PHENO_GROUP_CASE to the two groups you want "
            "to contrast, or narrow GROUP_DEFINITION so only those two match samples."
        )
    if requested_mode in {"temporal", "both"} and not can_run_temporal:
        raise ValueError(
            f"ANALYSIS_MODE='{requested_mode}' requires TEMPORAL_CORRELATION.enabled: true "
            "in the config."
        )
    if requested_mode == "auto" and not (can_run_pheno or can_run_temporal):
        raise ValueError(
            "ANALYSIS_MODE='auto' found nothing to run: the phenotype contrast needs exactly "
            f"two groups with samples (found {len(unique_groups)}: {unique_groups}), and "
            "TEMPORAL_CORRELATION is not enabled.\n"
            "To resolve: set PHENO_GROUP_REF/PHENO_GROUP_CASE, narrow GROUP_DEFINITION, or "
            "enable TEMPORAL_CORRELATION."
        )

    # Sample size guards — applied to all groups that will be used
    _check_group_sample_sizes(meta2["group"])

    # Connectivity-based selection may add features the abundance step did not
    # find, so the hard stop is applied after the two are combined.
    abundance_hard_stop = 0
    if not config.NULL_MODEL.enabled:
        logger.warning(
            "NULL_MODEL is disabled: networks and rewiring will be reported without any "
            "comparison to chance."
        )

    report_data: dict = {
        "tool_name": "PhenoRewire",
        "version": __version__,
        "config_path": rel(Path(config_path)) if config_path is not None else None,
        "null_model_enabled": bool(config.NULL_MODEL.enabled),
        "analysis_mode": requested_mode,
        "outdir": str(outdir),
        "n_samples_used": int(meta2.shape[0]),
        "n_features_input": int(len(df)),
        "n_features_after_preprocessing": int(X_corr.shape[0]),
        "log_transform": bool(config.LOG_TRANSFORM),
        "norm_method": str(config.NORM_METHOD),
        "triage_topology_weight": float(config.TRIAGE_TOPOLOGY_WEIGHT),
        "triage_selection_weight": float(config.TRIAGE_SELECTION_WEIGHT),
        "random_seed": int(config.RANDOM_SEED),
        "input_check": {"warnings": check.warnings},
        "deliverables": [],
        "phenotype": {"executed": False},
        "temporal": {"executed": False},
    }

    if can_run_pheno:
        if ref_group is None or case_group is None:
            logger.warning("PHENO_GROUP_REF/PHENO_GROUP_CASE not set. Falling back to alphabetical assignment.")
            ref_group, case_group = unique_groups[0], unique_groups[1]
        else:
            ref_group = str(ref_group)
            case_group = str(case_group)
            if set([ref_group, case_group]) != set(unique_groups):
                raise ValueError(
                    f"PHENO_GROUP_REF/CASE ({ref_group}, {case_group}) do not match groups in data: {unique_groups}"
                )

        group_to_bin = {ref_group: 0, case_group: 1}
        y = meta2["group"].map(group_to_bin).values.astype(int)
        report_data["phenotype"].update(
            {"executed": True, "ref_group": ref_group, "case_group": case_group}
        )
    else:
        # Only reachable in 'auto' mode with a temporal analysis to fall back on;
        # every other combination has already raised above.
        logger.debug(
            "Phenotype comparison not run: %d group(s) with samples after filtering (%s); "
            "continuing with the temporal analysis.",
            len(unique_groups),
            unique_groups,
        )
        y = np.array([], dtype=int)

    # -----------------------------
    # A) Phenotype comparison
    # -----------------------------
    if can_run_pheno:
        comp = layout.phenotype_comparison(outdir, ref_group, case_group)
        report_data["phenotype"]["folder"] = rel(comp.root)
        pheno_label = f"{ref_group} vs {case_group}"
        steps.start(f"Selecting metabolites: {pheno_label}")
        pheno_out = ensure_dir(comp.selection)
        pheno_all, selected, pheno_sel_summary = phenotype_selection(
            X=X_corr,
            y=y,
            feat_anno=feat_anno,
            outdir=pheno_out,
            n_permutations=config.N_PERMUTATIONS_PHENO,
            fdr_alpha=config.FDR_ALPHA,
            adaptive_fdr=config.ADAPTIVE_SELECTION_THRESHOLDS,
            min_selected_features=config.SELECTION_MIN_FEATURES,
            fdr_alpha_ceiling=config.SELECTION_FDR_ALPHA_CEILING,
            seed=config.RANDOM_SEED,
            min_features_hard_stop=abundance_hard_stop,
            permutation_stream=PERMUTATION_STREAM,
        )

        selected, selection_meta = _apply_connectivity_selection(
            config,
            pheno_all,
            selected,
            X_corr.loc[:, _subset_samples(meta2, config.META_SAMPLE_COL, meta2["group"] == ref_group)],
            X_corr.loc[:, _subset_samples(meta2, config.META_SAMPLE_COL, meta2["group"] == case_group)],
            pheno_out,
            q_col="q_pheno",
            label=f"Phenotype comparison {pheno_label}",
        )
        pheno_sel_summary["connectivity"] = selection_meta
        feats_sel = selected.index.astype(str).tolist()
        report_data["phenotype"]["n_selected_features"] = int(len(feats_sel))
        report_data["phenotype"]["selection"] = pheno_sel_summary
        steps.detail(_selection_line(len(feats_sel), selection_meta, "abundance"))
        if len(feats_sel) == 0:
            logger.warning("No phenotype-associated features selected. Skipping phenotype networks.")
            steps.skip(2)
        else:
            if not config.LOG_TRANSFORM:
                logger.warning(
                    "LOG_TRANSFORM is False: 'log2FC_case_vs_ref' contains raw "
                    "mean differences, not log2 fold-changes. Interpret with caution."
                )
            selected_pheno = _add_directional_columns(
                selected=selected,
                X_corr=X_corr,
                meta2=meta2,
                sample_col=config.META_SAMPLE_COL,
                ref_group=ref_group,
                case_group=case_group,
            )

            samp_ref = _subset_samples(meta2, config.META_SAMPLE_COL, meta2["group"] == ref_group)
            samp_case = _subset_samples(meta2, config.META_SAMPLE_COL, meta2["group"] == case_group)

            op_range = operating_range(len(samp_ref), len(samp_case))
            _log_operating_range(f"Phenotype comparison {pheno_label}", op_range)
            report_data["phenotype"]["operating_range"] = op_range

            report_ref, report_case, rewiring = _compare_states(
                config,
                comp,
                selected_pheno,
                X_corr.loc[feats_sel, samp_ref],
                X_corr.loc[feats_sel, samp_case],
                meta2.assign(state=meta2["group"]),
                kind="phenotype",
                label_a=ref_group,
                label_b=case_group,
                stop_labels=(ref_group, case_group),
                direction_col="direction",
                steps=steps,
            )
            report_ref["label"] = ref_group
            report_case["label"] = case_group
            report_data["phenotype"]["network_ref"] = report_ref
            report_data["phenotype"]["network_case"] = report_case
            report_data["phenotype"]["rewiring"] = rewiring
            report_data["deliverables"].extend(_key_files(comp, rewiring))

    # -----------------------------
    # B) Temporal comparison (optional)
    # -----------------------------
    if can_run_temporal:
        pheno_name = config.TEMPORAL_CORRELATION.phenotype
        time_col = config.TEMPORAL_CORRELATION.time_column

        if pheno_name is None or time_col is None:
            raise ValueError("TEMPORAL_CORRELATION.phenotype and time_column must be set when enabled.")

        meta_sub = meta2[meta2["group"] == pheno_name].copy()
        if meta_sub.empty:
            logger.warning("No samples found for phenotype '%s' in temporal analysis.", pheno_name)
            steps.skip(3)
        else:
            min_samp = int(getattr(config.TEMPORAL_CORRELATION, "min_samples", 3))
            if meta_sub.shape[0] < min_samp:
                logger.warning(
                    "Temporal analysis skipped for phenotype '%s': only %d samples (< min_samples=%d).",
                    pheno_name, meta_sub.shape[0], min_samp
                )
                steps.skip(3)
            else:
                if time_col not in meta_sub.columns:
                    raise ValueError(f"Time column '{time_col}' not found in metadata.")

                time_num = parse_time_numeric(meta_sub[time_col])
                if time_num.isna().any():
                    bad = meta_sub.loc[time_num.isna(), [config.META_SAMPLE_COL, time_col]].head(10)
                    raise ValueError(
                        "Could not parse some time values to numeric. Examples:\n"
                        f"{bad.to_string(index=False)}"
                    )

                t1, t2 = float(config.TIMEPOINT_T1), float(config.TIMEPOINT_T2)
                lbl_t1, lbl_t2 = timepoint_label(t1), timepoint_label(t2)
                comp_t = layout.temporal_comparison(outdir, pheno_name, lbl_t1, lbl_t2)
                report_data["temporal"]["folder"] = rel(comp_t.root)
                steps.start(f"Selecting metabolites: {pheno_name}, {lbl_t1} vs {lbl_t2}")

                samp_ids_sub = meta_sub[config.META_SAMPLE_COL].astype(str).tolist()
                X_sub = X_corr.loc[:, samp_ids_sub]

                temp_out = ensure_dir(comp_t.selection)
                temp_all, temp_selected, temp_sel_summary = temporal_selection(
                    X=X_sub,
                    time=time_num.values.astype(float),
                    feat_anno=feat_anno,
                    outdir=temp_out,
                    n_permutations=config.TEMPORAL_CORRELATION.n_permutations,
                    fdr_alpha=config.TEMPORAL_CORRELATION.fdr_alpha,
                    adaptive_fdr=config.ADAPTIVE_SELECTION_THRESHOLDS,
                    min_selected_features=config.SELECTION_MIN_FEATURES,
                    fdr_alpha_ceiling=config.SELECTION_FDR_ALPHA_CEILING,
                    seed=config.RANDOM_SEED,
                    min_features_hard_stop=abundance_hard_stop,
                    permutation_stream=PERMUTATION_STREAM,
                )

                temp_selected, selection_meta_t = _apply_connectivity_selection(
                    config,
                    temp_all,
                    temp_selected,
                    X_corr.loc[:, meta_sub.loc[(time_num == t1).values, config.META_SAMPLE_COL].astype(str)],
                    X_corr.loc[:, meta_sub.loc[(time_num == t2).values, config.META_SAMPLE_COL].astype(str)],
                    temp_out,
                    q_col="q_time",
                    label=f"Temporal comparison {lbl_t1} vs {lbl_t2}",
                )
                temp_sel_summary["connectivity"] = selection_meta_t

                temp_feats = temp_selected.index.astype(str).tolist()
                report_data["temporal"].update(
                    {
                        "executed": True,
                        "phenotype": pheno_name,
                        "n_selected_features": int(len(temp_feats)),
                        "selection": temp_sel_summary,
                    }
                )
                steps.detail(_selection_line(len(temp_feats), selection_meta_t, "association with time"))
                if len(temp_feats) == 0:
                    logger.warning("No temporal features selected. Skipping temporal rewiring networks.")
                    steps.skip(2)
                else:
                    meta_t1 = meta_sub[time_num == t1].copy()
                    meta_t2 = meta_sub[time_num == t2].copy()
                    other_tp = sorted(set(time_num.dropna().unique()) - {t1, t2})
                    if other_tp:
                        steps.detail(
                            "networks compare " + f"{lbl_t1} and {lbl_t2}; the other timepoints ("
                            + ", ".join(timepoint_label(t) for t in other_tp)
                            + ") were used for selection only"
                        )
                    if meta_t1.empty or meta_t2.empty:
                        raise ValueError(
                            f"Temporal rewiring requires samples at both timepoints. "
                            f"Found n(t1={t1})={meta_t1.shape[0]} and n(t2={t2})={meta_t2.shape[0]} "
                            f"for phenotype='{pheno_name}'."
                        )

                    samp_t1 = meta_t1[config.META_SAMPLE_COL].astype(str).tolist()
                    samp_t2 = meta_t2[config.META_SAMPLE_COL].astype(str).tolist()

                    X_t1 = X_corr.loc[temp_feats, samp_t1]
                    X_t2 = X_corr.loc[temp_feats, samp_t2]

                    op_range_t = operating_range(len(samp_t1), len(samp_t2))
                    _log_operating_range(f"Temporal comparison {lbl_t1} vs {lbl_t2}", op_range_t)
                    report_data["temporal"]["operating_range"] = op_range_t

                    temp_sel2 = temp_selected.copy()
                    mean_t1 = _safe_mean_rows(X_t1)
                    mean_t2 = _safe_mean_rows(X_t2)
                    temp_sel2["mean_t1"] = mean_t1.reindex(temp_feats).values
                    temp_sel2["mean_t2"] = mean_t2.reindex(temp_feats).values
                    temp_sel2["delta_t2_minus_t1"] = temp_sel2["mean_t2"].astype(float) - temp_sel2["mean_t1"].astype(float)

                    eps = 1e-12
                    temp_sel2["time_direction"] = np.where(
                        temp_sel2["delta_t2_minus_t1"].astype(float) > eps,
                        f"up_at_{lbl_t2}",
                        np.where(temp_sel2["delta_t2_minus_t1"].astype(float) < -eps, f"up_at_{lbl_t1}", "flat"),
                    )

                    samples_t = pd.concat([meta_t1.assign(state=lbl_t1), meta_t2.assign(state=lbl_t2)])
                    report_t1, report_t2, rewiring_t = _compare_states(
                        config,
                        comp_t,
                        temp_sel2,
                        X_t1,
                        X_t2,
                        samples_t,
                        kind="temporal",
                        label_a=lbl_t1,
                        label_b=lbl_t2,
                        stop_labels=(f"t={lbl_t1}", f"t={lbl_t2}"),
                        direction_col="time_direction",
                        steps=steps,
                    )
                    report_t1["label"] = f"{lbl_t1}"
                    report_t2["label"] = f"{lbl_t2}"
                    report_data["temporal"].update(
                        {
                            "timepoint_a": f"{lbl_t1}",
                            "timepoint_b": f"{lbl_t2}",
                            "network_t1": report_t1,
                            "network_t2": report_t2,
                            "rewiring": rewiring_t,
                        }
                    )
                    report_data["deliverables"].extend(_key_files(comp_t, rewiring_t))

    steps.start("Writing the report")
    report_data["deliverables"].append(str(outdir / layout.REPORT_MD))
    export_final_report(report_data, outdir)

    _closing_block(report_data, outdir, time.time() - started)


def _key_files(comp: layout.Comparison, rewiring: dict) -> list[str]:
    files = [comp.rewiring_vs_chance] if rewiring.get("null") else []
    files += [comp.rewired_metabolites, comp.selected_metabolites, comp.rewiring_network]
    return [str(p) for p in files]


def _closing_block(report: dict, outdir: Path, seconds: float) -> None:
    progress.info("")
    progress.info("Done in %s.", _elapsed(seconds))
    for key, a_key, b_key in (("phenotype", "ref_group", "case_group"), ("temporal", "timepoint_a", "timepoint_b")):
        section = report.get(key, {})
        if section.get("rewiring"):
            verdict = comparison_verdict(section, section.get(a_key), section.get(b_key))
            Steps.detail(f"{section.get(a_key)} vs {section.get(b_key)}: {verdict}", indent="  ")
    progress.info("")
    progress.info("  Start here : %s", rel(outdir / layout.REPORT_MD))
    for key in ("phenotype", "temporal"):
        folder = report.get(key, {}).get("folder")
        if folder and report[key].get("rewiring"):
            progress.info("  Key tables : %s/", folder)
    progress.info("  Run log    : %s", rel(layout.run_info(outdir).log))
