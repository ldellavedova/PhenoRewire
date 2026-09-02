from __future__ import annotations

from pathlib import Path
import logging
import shutil

import numpy as np
import pandas as pd

from .config import LeanConfig
from .utils import ensure_dir
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
from .temporal_selection import temporal_selection
from .triage import export_network_triage, export_rewiring_triage
from .reporting import export_final_report

logger = logging.getLogger(__name__)


# ============================================================
# Logging
# ============================================================

def _setup_logging(level: str = "INFO") -> None:
    lvl = getattr(logging, str(level).upper(), logging.INFO)
    logging.basicConfig(
        level=lvl,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )


# ============================================================
# Time helpers
# ============================================================

def _parse_time_numeric(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip().str.lower()
    s = s.str.replace("hours", "", regex=False).str.replace("hour", "", regex=False).str.replace("h", "", regex=False)
    s = s.str.strip()
    return pd.to_numeric(s, errors="coerce")


def _timepoint_label(t: float) -> str:
    """Filename-safe label for a timepoint.

    Integer-valued timepoints keep the plain "3h" form.  Fractional ones keep their
    value with the decimal point spelled out ("0.5" -> "0p5h"), so two distinct
    timepoints can never collapse onto the same directory name.
    """
    t = float(t)
    if t.is_integer():
        return f"{int(t)}h"
    return f"{t:g}".replace("-", "neg").replace(".", "p") + "h"


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
    community_summary.to_csv(outdir / "community_summary.csv", index=False)
    network_summary.to_csv(outdir / "network_summary.csv", index=False)
    export_network_triage(
        nodes_annotated,
        community_summary,
        network_summary,
        outdir,
        top_n=top_n,
        topology_weight=topology_weight,
        selection_weight=selection_weight,
    )

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

    priority_path = outdir / "triage" / "priority_features.csv"
    priority_df = pd.read_csv(priority_path) if priority_path.exists() else pd.DataFrame()
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
        "priority_features": priority_df.head(top_n).to_dict(orient="records"),
        "priority_columns": priority_df.columns.tolist(),
        "louvain_stability": louvain_stability_meta,
    }
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


def _copy_to_run_triage(paths: list[Path], triage_dir: Path) -> None:
    triage_dir.mkdir(parents=True, exist_ok=True)
    for path in paths:
        src = Path(path)
        if src.exists():
            shutil.copy2(src, triage_dir / src.name)


# ============================================================
# Main
# ============================================================

def run(config: LeanConfig) -> None:
    logger.info("Starting PhenoRewire...")

    df = pd.read_csv(config.DATA_MATRIX, low_memory=False)
    meta = pd.read_csv(config.METADATA, low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]
    meta.columns = [str(c).strip() for c in meta.columns]

    outdir = ensure_dir(Path(config.OUTDIR))
    run_triage_dir = ensure_dir(outdir / "triage")
    top_n = int(config.REPORT_TOP_N)
    requested_mode = str(config.ANALYSIS_MODE).lower()

    # -----------------------------
    # Preprocessing -> matrices
    # -----------------------------
    X_corr, meta2, feat_anno = prepare_matrices(
        df=df,
        meta=meta,
        feature_id_col=config.FEATURE_ID_COL,
        mz_col=config.MZ_COL,
        rt_col=config.RT_COL,
        name_col=config.NAME_COL,
        intensity_regex=config.INTENSITY_REGEX,
        meta_sample_col=config.META_SAMPLE_COL,
        outdir=outdir,
        min_nonzero_intensity=config.MIN_NONZERO_INTENSITY,
        min_presence_global=config.MIN_PRESENCE_GLOBAL,
        norm_method=config.NORM_METHOD,
        log_transform=config.LOG_TRANSFORM,
        group_definition=config.GROUP_DEFINITION,
    )

    meta2["group"].value_counts(dropna=False).rename_axis("group").to_frame("n").to_csv(outdir / "group_counts_used.csv")
    logger.info("Groups used:\n%s", meta2["group"].value_counts(dropna=False).to_string())

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
        logger.info(
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

    report_data: dict = {
        "tool_name": "PhenoRewire",
        "analysis_mode": requested_mode,
        "outdir": str(outdir),
        "n_samples_used": int(meta2.shape[0]),
        "n_features_after_preprocessing": int(X_corr.shape[0]),
        "triage_topology_weight": float(config.TRIAGE_TOPOLOGY_WEIGHT),
        "triage_selection_weight": float(config.TRIAGE_SELECTION_WEIGHT),
        "random_seed": int(config.RANDOM_SEED),
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

        pd.DataFrame([{"group_ref": ref_group, "group_case": case_group, "ref_bin": 0, "case_bin": 1}]).to_csv(
            outdir / "phenotype_direction.csv", index=False
        )
        report_data["phenotype"].update(
            {"executed": True, "ref_group": ref_group, "case_group": case_group}
        )
    else:
        # Only reachable in 'auto' mode with a temporal analysis to fall back on;
        # every other combination has already raised above.
        logger.info(
            "Phenotype comparison not run: %d group(s) with samples after filtering (%s); "
            "continuing with the temporal analysis.",
            len(unique_groups),
            unique_groups,
        )
        y = np.array([], dtype=int)

    # -----------------------------
    # A) Phenotype selection
    # -----------------------------
    if can_run_pheno:
        pheno_out = ensure_dir(outdir / "phenotype_selection")
        _, selected, pheno_sel_summary = phenotype_selection(
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
            min_features_hard_stop=config.MIN_FEATURES_HARD_STOP,
        )

        feats_sel = selected.index.astype(str).tolist()
        report_data["phenotype"]["n_selected_features"] = int(len(feats_sel))
        report_data["phenotype"]["selection"] = pheno_sel_summary
        if len(feats_sel) == 0:
            logger.warning("No phenotype-associated features selected. Skipping phenotype networks.")
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

            X_sel_ref = X_corr.loc[feats_sel, samp_ref]
            X_sel_case = X_corr.loc[feats_sel, samp_case]

            net_ref_out = outdir / f"network_{ref_group}"
            net_case_out = outdir / f"network_{case_group}"

            edges_ref, report_ref = _build_and_export_network(
                selected_pheno,
                X_sel_ref,
                net_ref_out,
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
            )
            _check_network_hard_stops(report_ref, ref_group, config.NETWORK_MIN_EDGES_HARD_STOP, config.NETWORK_MIN_NODES_HARD_STOP)
            edges_case, report_case = _build_and_export_network(
                selected_pheno,
                X_sel_case,
                net_case_out,
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
            )
            _check_network_hard_stops(report_case, case_group, config.NETWORK_MIN_EDGES_HARD_STOP, config.NETWORK_MIN_NODES_HARD_STOP)
            report_ref["label"] = ref_group
            report_case["label"] = case_group
            report_data["phenotype"]["network_ref"] = report_ref
            report_data["phenotype"]["network_case"] = report_case
            report_data["deliverables"].extend(
                [
                    str(net_ref_out / "triage" / "priority_features.csv"),
                    str(net_case_out / "triage" / "priority_features.csv"),
                ]
            )

            rew_pheno_out = ensure_dir(outdir / "rewiring_pheno")
            shared, only_ref, only_case, union = compute_rewiring(
                edges_ref,
                edges_case,
                label_a=ref_group,
                label_b=case_group,
                sign_switch_min_r=config.SIGN_SWITCH_MIN_R,
            )
            rewiring_summary, rewiring_nodes = summarize_rewiring(
                union,
                label_a=ref_group,
                label_b=case_group,
            )
            shared.to_csv(rew_pheno_out / "edges_shared.csv", index=False)
            only_ref.to_csv(rew_pheno_out / f"edges_{ref_group}_only.csv", index=False)
            only_case.to_csv(rew_pheno_out / f"edges_{case_group}_only.csv", index=False)
            union.to_csv(rew_pheno_out / "edges_union_with_stats.csv", index=False)
            rewiring_summary.to_csv(rew_pheno_out / "rewiring_summary.csv", index=False)
            rewiring_nodes.to_csv(rew_pheno_out / "rewiring_node_summary.csv", index=False)
            meta2.to_csv(rew_pheno_out / "metadata_used.csv", index=False)
            selected_pheno.to_csv(pheno_out / "selected_pheno_features_with_direction.csv")
            export_rewiring_triage(rewiring_summary, rewiring_nodes, rew_pheno_out, top_n=top_n)
            _pheno_graphml = outdir / "triage" / "rewiring_network_phenotype.graphml"
            export_rewiring_network(selected_pheno, union, rewiring_nodes, _pheno_graphml)
            logger.info("Rewiring network saved: %s", _pheno_graphml)
            report_data["phenotype"]["rewiring"] = {
                "summary": rewiring_summary.iloc[0].to_dict() if not rewiring_summary.empty else {},
                "priority_nodes": rewiring_nodes.head(top_n).to_dict(orient="records"),
            }
            report_data["deliverables"].extend(
                [
                    str(pheno_out / "selected_pheno_features_with_direction.csv"),
                    str(rew_pheno_out / "triage" / "priority_rewired_nodes.csv"),
                    str(_pheno_graphml),
                ]
            )
            _copy_to_run_triage(
                [
                    pheno_out / "selected_pheno_features_with_direction.csv",
                    rew_pheno_out / "triage" / "priority_rewired_nodes.csv",
                    rew_pheno_out / "triage" / "rewiring_summary.csv",
                ],
                run_triage_dir,
            )

    # -----------------------------
    # B) Temporal selection + rewiring (optional)
    # -----------------------------
    if can_run_temporal:
        logger.info("Temporal correlation enabled.")
        pheno_name = config.TEMPORAL_CORRELATION.phenotype
        time_col = config.TEMPORAL_CORRELATION.time_column

        if pheno_name is None or time_col is None:
            raise ValueError("TEMPORAL_CORRELATION.phenotype and time_column must be set when enabled.")

        meta_sub = meta2[meta2["group"] == pheno_name].copy()
        if meta_sub.empty:
            logger.warning("No samples found for phenotype '%s' in temporal analysis.", pheno_name)
        else:
            min_samp = int(getattr(config.TEMPORAL_CORRELATION, "min_samples", 3))
            if meta_sub.shape[0] < min_samp:
                logger.warning(
                    "Temporal analysis skipped for phenotype '%s': only %d samples (< min_samples=%d).",
                    pheno_name, meta_sub.shape[0], min_samp
                )
            else:
                if time_col not in meta_sub.columns:
                    raise ValueError(f"Time column '{time_col}' not found in metadata.")

                time_num = _parse_time_numeric(meta_sub[time_col])
                if time_num.isna().any():
                    bad = meta_sub.loc[time_num.isna(), [config.META_SAMPLE_COL, time_col]].head(10)
                    raise ValueError(
                        "Could not parse some time values to numeric. Examples:\n"
                        f"{bad.to_string(index=False)}"
                    )

                samp_ids_sub = meta_sub[config.META_SAMPLE_COL].astype(str).tolist()
                X_sub = X_corr.loc[:, samp_ids_sub]

                temp_out = ensure_dir(outdir / "temporal_selection")
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
                    min_features_hard_stop=config.MIN_FEATURES_HARD_STOP,
                )

                temp_feats = temp_selected.index.astype(str).tolist()
                report_data["temporal"].update(
                    {
                        "executed": True,
                        "phenotype": pheno_name,
                        "n_selected_features": int(len(temp_feats)),
                        "selection": temp_sel_summary,
                    }
                )
                if len(temp_feats) == 0:
                    logger.warning("No temporal features selected. Skipping temporal rewiring networks.")
                else:
                    t1 = float(getattr(config, "TIMEPOINT_T1", 3))
                    lbl_t1 = _timepoint_label(t1)
                    t2 = float(getattr(config, "TIMEPOINT_T2", 6))
                    lbl_t2 = _timepoint_label(t2)

                    meta_t1 = meta_sub[time_num == t1].copy()
                    meta_t2 = meta_sub[time_num == t2].copy()
                    other_tp = sorted(set(time_num.dropna().unique()) - {t1, t2})
                    if other_tp:
                        logger.info(
                            "Temporal analysis: using T1=%.4g and T2=%.4g. "
                            "Ignoring other timepoints present in '%s': %s",
                            t1, t2, pheno_name, other_tp,
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

                    net_t1_out = outdir / f"network_{pheno_name}_{lbl_t1}"
                    net_t2_out = outdir / f"network_{pheno_name}_{lbl_t2}"

                    edges_t1, report_t1 = _build_and_export_network(
                        temp_sel2,
                        X_t1,
                        net_t1_out,
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
                    )
                    _check_network_hard_stops(report_t1, f"t={lbl_t1}", config.NETWORK_MIN_EDGES_HARD_STOP, config.NETWORK_MIN_NODES_HARD_STOP)
                    edges_t2, report_t2 = _build_and_export_network(
                        temp_sel2,
                        X_t2,
                        net_t2_out,
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
                    )
                    _check_network_hard_stops(report_t2, f"t={lbl_t2}", config.NETWORK_MIN_EDGES_HARD_STOP, config.NETWORK_MIN_NODES_HARD_STOP)
                    report_t1["label"] = f"{lbl_t1}"
                    report_t2["label"] = f"{lbl_t2}"
                    report_data["temporal"].update(
                        {
                            "timepoint_a": f"{lbl_t1}",
                            "timepoint_b": f"{lbl_t2}",
                            "network_t1": report_t1,
                            "network_t2": report_t2,
                        }
                    )
                    report_data["deliverables"].extend(
                        [
                            str(net_t1_out / "triage" / "priority_features.csv"),
                            str(net_t2_out / "triage" / "priority_features.csv"),
                        ]
                    )

                    rew_time_out = ensure_dir(outdir / "rewiring_time")
                    shared_t, only_t1, only_t2, union_t = compute_rewiring(
                        edges_t1,
                        edges_t2,
                        label_a=f"{lbl_t1}",
                        label_b=f"{lbl_t2}",
                        sign_switch_min_r=config.SIGN_SWITCH_MIN_R,
                    )
                    rewiring_summary_t, rewiring_nodes_t = summarize_rewiring(
                        union_t,
                        label_a=f"{lbl_t1}",
                        label_b=f"{lbl_t2}",
                    )
                    shared_t.to_csv(rew_time_out / "edges_shared.csv", index=False)
                    only_t1.to_csv(rew_time_out / f"edges_{lbl_t1}_only.csv", index=False)
                    only_t2.to_csv(rew_time_out / f"edges_{lbl_t2}_only.csv", index=False)
                    union_t.to_csv(rew_time_out / "edges_union_with_stats.csv", index=False)
                    rewiring_summary_t.to_csv(rew_time_out / "rewiring_summary.csv", index=False)
                    rewiring_nodes_t.to_csv(rew_time_out / "rewiring_node_summary.csv", index=False)

                    meta_sub.to_csv(rew_time_out / "metadata_used.csv", index=False)
                    temp_sel2.to_csv(temp_out / "selected_temporal_features_with_direction.csv")
                    export_rewiring_triage(rewiring_summary_t, rewiring_nodes_t, rew_time_out, top_n=top_n)
                    _temp_graphml = outdir / "triage" / "rewiring_network_temporal.graphml"
                    export_rewiring_network(temp_sel2, union_t, rewiring_nodes_t, _temp_graphml)
                    logger.info("Rewiring network saved: %s", _temp_graphml)
                    report_data["temporal"]["rewiring"] = {
                        "summary": rewiring_summary_t.iloc[0].to_dict() if not rewiring_summary_t.empty else {},
                        "priority_nodes": rewiring_nodes_t.head(top_n).to_dict(orient="records"),
                    }
                    report_data["deliverables"].extend(
                        [
                            str(temp_out / "selected_temporal_features_with_direction.csv"),
                            str(rew_time_out / "triage" / "priority_rewired_nodes.csv"),
                            str(_temp_graphml),
                        ]
                    )
                    _copy_to_run_triage(
                        [
                            temp_out / "selected_temporal_features_with_direction.csv",
                            rew_time_out / "triage" / "priority_rewired_nodes.csv",
                            rew_time_out / "triage" / "rewiring_summary.csv",
                        ],
                        run_triage_dir,
                    )

    report_data["deliverables"].append(str(outdir / "report" / "final_report.md"))
    export_final_report(report_data, outdir)
    _copy_to_run_triage(
        [
            outdir / "report" / "final_report.md",
            outdir / "report" / "final_report.json",
        ],
        run_triage_dir,
    )

    logger.info("Pipeline completed. OUTDIR = %s", outdir)
