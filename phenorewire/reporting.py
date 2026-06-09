# ── PhenoRewire · Step 5 — Reporting ─────────────────────────────────────────
# Generates final_report.md and final_report.json summarising the full run,
# including selected features, network statistics, and rewiring summary.
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_dir, save_json


def _relativize_path(path_str: str) -> str:
    """Return path relative to cwd when possible; always use forward slashes."""
    try:
        return str(Path(str(path_str)).relative_to(Path.cwd())).replace("\\", "/")
    except ValueError:
        return str(path_str).replace("\\", "/")


def _fmt_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _top_records(df: pd.DataFrame | None, cols: list[str], top_n: int) -> list[dict[str, Any]]:
    if df is None or df.empty:
        return []
    keep = [col for col in cols if col in df.columns]
    return df.loc[:, keep].head(top_n).to_dict(orient="records")


def _bullet_lines(records: list[dict[str, Any]], primary_col: str) -> list[str]:
    lines: list[str] = []
    for row in records:
        lead = row.get(primary_col, "NA")
        rest = [f"{k}={_fmt_value(v)}" for k, v in row.items() if k != primary_col]
        if rest:
            lines.append(f"- {lead}: " + ", ".join(rest))
        else:
            lines.append(f"- {lead}")
    return lines


def _collect_warnings(report: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    phenotype = report.get("phenotype", {})
    temporal = report.get("temporal", {})

    if phenotype.get("executed") and int(phenotype.get("n_selected_features", 0)) == 0:
        warnings.append("Phenotype comparison selected no significant features under the current thresholds.")
    if temporal.get("executed") and int(temporal.get("n_selected_features", 0)) == 0:
        warnings.append("Temporal analysis selected no significant features under the current thresholds.")

    # Adaptive FDR warnings — prominent because reviewers will ask about this
    for section_name, section in [("Phenotype", phenotype), ("Temporal", temporal)]:
        sel = section.get("selection", {})
        if sel and sel.get("adaptive_fdr_used"):
            base_a = sel.get("fdr_alpha", "?")
            applied_a = sel.get("applied_fdr_alpha", "?")
            warnings.append(
                f"⚠ ADAPTIVE FDR RELAXATION ({section_name}): the configured FDR alpha "
                f"({base_a}) yielded too few features; the threshold was relaxed to "
                f"{applied_a}. Results from the relaxed threshold are flagged but should be "
                "interpreted with caution. Consider increasing sample size or explicitly "
                "setting ADAPTIVE_SELECTION_THRESHOLDS: false if you want strict control."
            )

    for section in [phenotype, temporal]:
        for key in ["network_ref", "network_case", "network_t1", "network_t2"]:
            net = section.get(key, {})
            if not net:
                continue
            if int(net.get("n_edges", 0)) == 0:
                warnings.append(
                    f"Network '{net.get('label', 'NA')}' has no significant edges after filtering. "
                    "Consider relaxing CORR_ABS_THRESHOLD, raising CORR_FDR_ALPHA, or adding more samples."
                )
            elif float(net.get("modularity", 0.0)) < 0.2:
                warnings.append(
                    f"Network '{net.get('label', 'NA')}' has low modularity ({net.get('modularity', 0):.3f}); "
                    "community structure may be weak or driven by a single dense cluster."
                )
            # Louvain stability warning
            louvain = net.get("louvain_stability", {})
            if louvain and not louvain.get("stable", True):
                warnings.append(
                    f"Community assignments for network '{net.get('label', 'NA')}' show low stability "
                    f"(mean AMI = {louvain.get('ami_mean', '?')} over {louvain.get('n_runs', '?')} runs). "
                    "Module memberships are not reliable; consider adjusting CORR_ABS_THRESHOLD."
                )

    # Sign-switch reporting
    for section_name, section in [("phenotype", phenotype), ("temporal", temporal)]:
        rew = section.get("rewiring", {})
        summary = rew.get("summary", {})
        if summary:
            n_switch = int(summary.get("n_sign_switch_edges", 0))
            n_shared = int(summary.get("n_shared_edges", 1))
            pct = 100.0 * n_switch / max(1, n_shared)
            if pct > 20:
                warnings.append(
                    f"⚠ High sign-switch rate ({section_name}): {n_switch} sign-switch edges "
                    f"out of {n_shared} shared edges ({pct:.1f}%). This may reflect marginal "
                    "correlations near zero. Check SIGN_SWITCH_MIN_R setting."
                )

    return warnings


def _collect_takeaways(report: dict[str, Any]) -> list[str]:
    notes: list[str] = []
    phenotype = report.get("phenotype", {})
    temporal = report.get("temporal", {})

    if phenotype.get("executed"):
        n_sel = int(phenotype.get("n_selected_features", 0))
        notes.append(f"Phenotype comparison retained {n_sel} phenotype-associated features.")
        rew = phenotype.get("rewiring", {}).get("summary", {})
        if rew:
            only_counts = [int(v) for k, v in rew.items() if k.endswith("_only_edges")]
            notes.append(
                "Phenotype rewiring summary: "
                f"shared={rew.get('n_shared_edges', 0)}, "
                f"condition_specific={sum(only_counts)}, "
                f"sign_switch={rew.get('n_sign_switch_edges', 0)}."
            )

    if temporal.get("executed"):
        n_sel = int(temporal.get("n_selected_features", 0))
        notes.append(f"Temporal analysis retained {n_sel} time-associated features.")
        rew = temporal.get("rewiring", {}).get("summary", {})
        if rew:
            notes.append(
                "Temporal rewiring summary: "
                f"shared={rew.get('n_shared_edges', 0)}, "
                f"timepoint_specific={sum(int(v) for k, v in rew.items() if k.endswith('_only_edges'))}, "
                f"sign_switch={rew.get('n_sign_switch_edges', 0)}."
            )
    return notes


def export_final_report(report: dict[str, Any], outdir: Path) -> None:
    outdir = ensure_dir(Path(outdir))
    report_dir = ensure_dir(outdir / "report")

    # Portable copy: relativize paths before writing user-facing reports
    report_portable = dict(report)
    report_portable["outdir"] = _relativize_path(str(report.get("outdir", "")))
    report_portable["deliverables"] = [
        _relativize_path(str(p)) for p in report.get("deliverables", [])
    ]
    save_json(report_portable, report_dir / "final_report.json")

    lines: list[str] = []
    title = report.get("tool_name", "PhenoRewire")
    lines.append(f"# {title} Final Report")
    lines.append("")
    lines.append("## Run Summary")
    lines.append(f"- analysis_mode: {report_portable.get('analysis_mode', 'NA')}")
    lines.append(f"- outdir: {report_portable.get('outdir', 'NA')}")
    lines.append(f"- n_samples_used: {report.get('n_samples_used', 'NA')}")
    lines.append(f"- n_features_after_preprocessing: {report.get('n_features_after_preprocessing', 'NA')}")
    lines.append("")
    lines.append("## Executive Summary")
    for line in _collect_takeaways(report):
        lines.append(f"- {line}")
    warning_lines = _collect_warnings(report)
    if not _collect_takeaways(report):
        lines.append("- Run completed successfully.")
    lines.append("")

    if warning_lines:
        lines.append("## Warnings")
        for line in warning_lines:
            lines.append(f"- {line}")
        lines.append("")

    phenotype = report.get("phenotype", {})
    if phenotype.get("executed"):
        lines.append("## Phenotype Comparison")
        lines.append(f"- groups: {phenotype.get('ref_group')} vs {phenotype.get('case_group')}")
        lines.append(f"- selected_features: {phenotype.get('n_selected_features', 0)}")
        if phenotype.get('selection'):
            sel = phenotype.get('selection', {})
            lines.append(f"- applied_selection_fdr_alpha: {_fmt_value(sel.get('applied_fdr_alpha', 'NA'))}")
            lines.append(f"- adaptive_selection_used: {sel.get('adaptive_fdr_used', False)}")
            lines.append(f"- selected_at_base_alpha: {sel.get('selected_at_base_alpha', 'NA')}")
        lines.append("")
        for key in ["network_ref", "network_case"]:
            net = phenotype.get(key, {})
            if not net:
                continue
            lines.append(f"### Network {net.get('label', 'NA')}")
            lines.append(f"- n_nodes: {net.get('n_nodes', 0)}")
            lines.append(f"- n_edges: {net.get('n_edges', 0)}")
            lines.append(f"- modularity: {_fmt_value(net.get('modularity', 0.0))}")
            lines.append(f"- applied_abs_r_thr: {_fmt_value(net.get('applied_abs_r_thr', 'NA'))}")
            lines.append(f"- applied_fdr_alpha: {_fmt_value(net.get('applied_fdr_alpha', 'NA'))}")
            lines.append(f"- adaptive_thresholds_used: {net.get('adaptive_thresholds_used', False)}")
            top_features = _top_records(pd.DataFrame(net.get("priority_features", [])), list(net.get("priority_columns", [])), 5)
            if top_features:
                lines.append("- top priority features:")
                lines.extend(_bullet_lines(top_features, "feature_id"))
            lines.append("")

        rew = phenotype.get("rewiring", {})
        if rew:
            lines.append("### Rewiring")
            for k, v in rew.get("summary", {}).items():
                lines.append(f"- {k}: {_fmt_value(v)}")
            top_nodes = rew.get("priority_nodes", [])
            if top_nodes:
                lines.append("- top rewired nodes:")
                lines.extend(_bullet_lines(top_nodes, "feature_id"))
            lines.append("")

    temporal = report.get("temporal", {})
    if temporal.get("executed"):
        lines.append("## Temporal Rewiring")
        lines.append(f"- phenotype: {temporal.get('phenotype', 'NA')}")
        lines.append(f"- timepoints: {temporal.get('timepoint_a', 'NA')} vs {temporal.get('timepoint_b', 'NA')}")
        lines.append(f"- selected_features: {temporal.get('n_selected_features', 0)}")
        if temporal.get('selection'):
            sel = temporal.get('selection', {})
            lines.append(f"- applied_selection_fdr_alpha: {_fmt_value(sel.get('applied_fdr_alpha', 'NA'))}")
            lines.append(f"- adaptive_selection_used: {sel.get('adaptive_fdr_used', False)}")
            lines.append(f"- selected_at_base_alpha: {sel.get('selected_at_base_alpha', 'NA')}")
        lines.append("")
        for key in ["network_t1", "network_t2"]:
            net = temporal.get(key, {})
            if not net:
                continue
            lines.append(f"### Network {net.get('label', 'NA')}")
            lines.append(f"- n_nodes: {net.get('n_nodes', 0)}")
            lines.append(f"- n_edges: {net.get('n_edges', 0)}")
            lines.append(f"- modularity: {_fmt_value(net.get('modularity', 0.0))}")
            lines.append(f"- applied_abs_r_thr: {_fmt_value(net.get('applied_abs_r_thr', 'NA'))}")
            lines.append(f"- applied_fdr_alpha: {_fmt_value(net.get('applied_fdr_alpha', 'NA'))}")
            lines.append(f"- adaptive_thresholds_used: {net.get('adaptive_thresholds_used', False)}")
            top_features = _top_records(pd.DataFrame(net.get("priority_features", [])), list(net.get("priority_columns", [])), 5)
            if top_features:
                lines.append("- top priority features:")
                lines.extend(_bullet_lines(top_features, "feature_id"))
            lines.append("")

        rew = temporal.get("rewiring", {})
        if rew:
            lines.append("### Rewiring")
            for k, v in rew.get("summary", {}).items():
                lines.append(f"- {k}: {_fmt_value(v)}")
            top_nodes = rew.get("priority_nodes", [])
            if top_nodes:
                lines.append("- top rewired nodes:")
                lines.extend(_bullet_lines(top_nodes, "feature_id"))
            lines.append("")

    _mode = str(report.get("analysis_mode", "phenotype")).lower()
    if _mode == "temporal":
        _graphml_name = "triage/rewiring_network_temporal.graphml"
    elif _mode == "both":
        _graphml_name = "triage/rewiring_network_phenotype.graphml (phenotype) or triage/rewiring_network_temporal.graphml (temporal)"
    else:
        _graphml_name = "triage/rewiring_network_phenotype.graphml"

    lines.append("## Start Here")
    lines.append("")
    lines.append("Open these four files in order:")
    lines.append("")
    lines.append("| # | File | What it tells you |")
    lines.append("| --- | --- | --- |")
    lines.append(f"| 1 | `{_graphml_name}` | Full rewiring network — open in Cytoscape; colour edges by `state`, size nodes by `rewiring_score` |")
    lines.append("| 2 | `triage/priority_rewired_nodes.csv` | Ranked list of the most rewired metabolites; start with the top 10–15 rows |")
    lines.append("| 3 | `triage/rewiring_summary.csv` | One-line global summary: how many edges are shared, state-specific, or sign-switched |")
    lines.append("| 4 | `triage/final_report.md` | This file — narrative summary of the full run |")
    lines.append("")
    lines.append("Everything else under `network_*/`, `phenotype_selection/`, and `rewiring_*/` is supporting detail.")
    lines.append("")
    t_w = float(report.get("triage_topology_weight", 0.7))
    s_w = float(report.get("triage_selection_weight", 0.3))
    rng_seed = report.get("random_seed", "NA")
    lines.append("## Scoring Formulas")
    lines.append("")
    lines.append("**triage_score** (per-network node ranking):")
    lines.append("```")
    lines.append(f"triage_score = {t_w} * topology_priority_score + {s_w} * selection_support_score")
    lines.append("topology_priority_score = mean percentile rank of [degree, betweenness, eigenvector, weighted_degree]")
    lines.append("selection_support_score = 1 - q_value  (higher = stronger phenotype/time association)")
    lines.append(f"  Configured: TRIAGE_TOPOLOGY_WEIGHT={t_w}, TRIAGE_SELECTION_WEIGHT={s_w}")
    lines.append("```")
    lines.append("")
    lines.append("**rewiring_score** (per-node rewiring ranking, absolute):")
    lines.append("```")
    lines.append("rewiring_score = state_a_only_edges + state_b_only_edges + sign_switch_edges")
    lines.append("```")
    lines.append("")
    lines.append("**rewiring_fraction** (per-node, normalised 0–1):")
    lines.append("```")
    lines.append("rewiring_fraction = (state_a_only + state_b_only) / total_union_degree")
    lines.append("Nodes with high rewiring_fraction rewire a large proportion of their connections,")
    lines.append("regardless of how many connections they have in total.")
    lines.append("```")
    lines.append("")
    lines.append("**Note — rewiring_score vs rewiring_fraction:**")
    lines.append("rewiring_score counts sign-switch edges as 'rewired'; rewiring_fraction does not.")
    lines.append("Use rewiring_score to rank 'most rewired' metabolites (absolute change including reversals).")
    lines.append("Use rewiring_fraction to compare metabolites with different degrees (normalised proportion).")
    lines.append("Sign-switch edges are listed separately in `edges_union_with_stats.csv` (state='sign_switch').")
    lines.append("")
    lines.append("**rewiring_proportion** (global summary):")
    lines.append("```")
    lines.append("rewiring_proportion = (n_a_only_edges + n_b_only_edges) / total_edges_in_union")
    lines.append("sign_switch_proportion = n_sign_switch_edges / n_shared_edges")
    lines.append("```")
    lines.append("")
    lines.append("## Reproducibility")
    lines.append("")
    lines.append("Every run should be fully reproducible from this block:")
    lines.append("```")
    lines.append(f"random_seed:              {rng_seed}")
    lines.append(f"triage_topology_weight:   {t_w}")
    lines.append(f"triage_selection_weight:  {s_w}")
    lines.append(f"analysis_mode:            {report.get('analysis_mode', 'NA')}")
    # Selection thresholds
    for section_name, section in [("phenotype", report.get("phenotype", {})), ("temporal", report.get("temporal", {}))]:
        sel = section.get("selection", {})
        if sel:
            lines.append(f"{section_name}_fdr_alpha:           {sel.get('fdr_alpha', 'NA')}")
            lines.append(f"{section_name}_applied_fdr_alpha:   {sel.get('applied_fdr_alpha', 'NA')}")
            lines.append(f"{section_name}_adaptive_fdr_used:   {sel.get('adaptive_fdr_used', False)}")
            lines.append(f"{section_name}_n_permutations:      {sel.get('n_permutations', 'NA')}")
    # Network thresholds
    for section_name, section in [("phenotype", report.get("phenotype", {})), ("temporal", report.get("temporal", {}))]:
        for net_key in ["network_ref", "network_case", "network_t1", "network_t2"]:
            net = section.get(net_key, {})
            if net:
                lbl = net.get("label", net_key)
                lines.append(f"network_{lbl}_applied_abs_r:  {net.get('applied_abs_r_thr', 'NA')}")
                lines.append(f"network_{lbl}_applied_fdr_q:  {net.get('applied_fdr_alpha', 'NA')}")
    lines.append("```")
    lines.append("")
    lines.append("## Next Steps")
    lines.append("")
    lines.append("1. **Targeted validation**: consider functional assays or literature lookup for the top 10 priority features.")
    lines.append("2. **Sign-switch edges**: metabolites connected by a sign-switch edge may indicate condition-specific metabolic competition; these are strong candidates for follow-up.")
    lines.append(f"3. **Cytoscape**: load `{_graphml_name}`, colour edges by `state`, size nodes by `rewiring_score`.")
    lines.append("4. **Interpretation caveat**: PhenoRewire is a hypothesis-generating tool; results are correlative, not causal.")
    lines.append("")
    lines.append("## All Output Files")
    for path in report_portable.get("deliverables", []):
        lines.append(f"- {path}")
    lines.append("")

    (report_dir / "final_report.md").write_text("\n".join(lines), encoding="utf-8")
