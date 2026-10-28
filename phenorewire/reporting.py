# ── PhenoRewire · Step 5 — Reporting ─────────────────────────────────────────
# Writes REPORT.md (what the run found, in words, and which file to open first)
# and report.json (everything the report says, machine-readable) at the top of
# the output folder.
from __future__ import annotations

from pathlib import Path
from typing import Any

from . import layout
from .console import fmt_p
from .utils import NORM_LABELS, ensure_dir, save_json


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


def _num(value: Any, digits: int = 2) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "NA"


def _rewiring_graphml_produced(report: dict[str, Any]) -> bool:
    """True when this run actually wrote a rewiring network, not just planned one."""
    return any(
        str(path).endswith(".graphml")
        for path in report.get("deliverables", [])
    )


def _sections(report: dict[str, Any]) -> list[tuple[str, dict[str, Any], str, str]]:
    """(name, section, state A, state B) for every analysis that ran."""
    out = []
    phenotype = report.get("phenotype", {})
    if phenotype.get("executed"):
        out.append(("Phenotype", phenotype, str(phenotype.get("ref_group")), str(phenotype.get("case_group"))))
    temporal = report.get("temporal", {})
    if temporal.get("executed"):
        out.append(("Temporal", temporal, str(temporal.get("timepoint_a")), str(temporal.get("timepoint_b"))))
    return out


def _title(name: str, section: dict[str, Any], a: str, b: str) -> str:
    if name == "Temporal":
        return f"{section.get('phenotype', '')}: {a} vs {b}"
    return f"{a} vs {b}"


def comparison_verdict(section: dict[str, Any], label_a: Any, label_b: Any) -> str:
    """The result of one comparison in one sentence."""
    rew = section.get("rewiring") or {}
    null = rew.get("null") or {}
    if null:
        obs = null.get("observed_state_specific_edges")
        expected = _num(null.get("null_mean"), 1)
        p = null.get("p_empirical")
        if p is not None and float(p) <= 0.05:
            n_nodes = int(null.get("n_nodes_q05", 0))
            per_node = (
                f"{n_nodes} metabolite(s) rewire more than chance (q <= 0.05)."
                if n_nodes
                else "no single metabolite does after correction for multiple testing."
            )
            return (
                f"the networks differ more than chance ({obs} state-specific correlations vs {expected} "
                f"expected, p = {fmt_p(p)}); {per_node}"
            )
        return (
            f"the difference between the networks is not distinguishable from chance ({obs} "
            f"state-specific correlations vs {expected} expected, p = {fmt_p(p)})."
        )
    summary = rew.get("summary") or {}
    if summary:
        n = summary.get("total_state_specific_edges", "NA")
        return f"{n} state-specific correlations; not compared to chance (NULL_MODEL is disabled)."
    return "no networks were compared."


def _collect_warnings(report: dict[str, Any]) -> list[str]:
    warnings: list[str] = list((report.get("input_check") or {}).get("warnings", []))
    phenotype = report.get("phenotype", {})
    temporal = report.get("temporal", {})

    if not phenotype.get("executed") and not temporal.get("executed"):
        warnings.append(
            "**No analysis was executed.** Neither the phenotype comparison nor the temporal "
            "analysis ran, so this report contains no results. Check ANALYSIS_MODE, "
            "GROUP_DEFINITION, and TEMPORAL_CORRELATION in your config."
        )
    elif not _rewiring_graphml_produced(report):
        warnings.append(
            "**No rewiring network was produced.** Feature selection ran but did not yield "
            "enough features to build and compare networks."
        )

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
                f"**Adaptive FDR relaxation ({section_name}):** the configured FDR alpha "
                f"({base_a}) yielded too few features; the threshold was relaxed to "
                f"{applied_a}. Results from the relaxed threshold should be "
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

    for section_name, section in [("Phenotype", phenotype), ("Temporal", temporal)]:
        op = section.get("operating_range") or {}
        if op.get("status") == "below":
            warnings.append(
                f"**Below the operating range ({section_name}):** {op.get('min_n_per_state')} samples in the "
                f"smaller state, against a minimum of {op.get('absolute_min_n')}. At this size networks "
                "are dominated by near-perfect rank correlations; do not interpret rewiring unless the "
                "rewiring null is clearly significant, and treat it as preliminary even then."
            )
        elif op.get("status") == "reduced":
            warnings.append(
                f"{section_name}: {op.get('min_n_per_state')} samples in the smaller state "
                f"(recommended >= {op.get('recommended_min_n')}). Recall is low at this size; an "
                "absence of rewiring is not evidence of no rewiring."
            )
        null = (section.get("rewiring") or {}).get("null") or {}
        p_val = null.get("p_empirical")
        if p_val is not None and float(p_val) > 0.05:
            warnings.append(
                f"**Rewiring not distinguishable from chance ({section_name}):** "
                f"{null.get('observed_state_specific_edges')} state-specific edges against "
                f"{_num(null.get('null_mean'), 1)} expected from an arbitrary split (p = {fmt_p(p_val)}). "
                "The per-metabolite ranking describes these data but should not be read as "
                "a difference between the states."
            )

    if report.get("null_model_enabled") is False:
        warnings.append(
            "NULL_MODEL was disabled for this run: no network or rewiring result has been compared "
            "to chance."
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
                    f"High sign-switch rate ({section_name}): {n_switch} sign-switch edges "
                    f"out of {n_shared} shared edges ({pct:.1f}%). This may reflect marginal "
                    "correlations near zero. Check SIGN_SWITCH_MIN_R setting."
                )

    return warnings


def _selection_sentence(section: dict[str, Any], what: str) -> str:
    sel = section.get("selection") or {}
    conn = sel.get("connectivity") or {}
    n = section.get("n_selected_features", 0)
    alpha = sel.get("applied_fdr_alpha", "NA")
    adaptive = " (relaxed adaptively)" if sel.get("adaptive_fdr_used") else ""
    if conn.get("connectivity_run"):
        return (
            f"{n} metabolites: {conn.get('n_abundance_only', 0)} by {what} only, "
            f"{conn.get('n_connectivity_only', 0)} only because their co-variation partners changed "
            f"(invisible to a differential-abundance analysis), {conn.get('n_both', 0)} by both; "
            f"FDR {alpha}{adaptive}."
        )
    reason = " Connectivity-based selection did not run: fewer than " \
        f"{conn.get('required_n', 14)} samples per state." if conn else ""
    return f"{n} metabolites by {what}, FDR {alpha}{adaptive}.{reason}"


def _network_rows(section: dict[str, Any], keys: tuple[str, str]) -> list[str]:
    rows = [
        "| Network | Samples | Metabolites | Correlations | Thresholds | Correlations expected by chance | Modularity beyond chance (ΔQ) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for key in keys:
        net = section.get(key) or {}
        if not net:
            continue
        null = net.get("null_model") or {}
        expected = f"{_num(null.get('null_mean'), 1)} (p = {fmt_p(null.get('p_empirical'))})" if null else "not tested"
        dq = _num(null.get("delta_modularity")) if null.get("delta_modularity") is not None else "NA"
        rows.append(
            f"| {net.get('label', key)} | {net.get('n_samples', 'NA')} | {net.get('n_nodes', 0)} | "
            f"{net.get('n_edges', 0)} | \\|r\\| ≥ {_num(net.get('applied_abs_r_thr'))}, "
            f"q ≤ {_num(net.get('applied_fdr_alpha'))} | {expected} | {dq} |"
        )
    return rows


def _top_rewired_rows(section: dict[str, Any]) -> list[str]:
    nodes = (section.get("rewiring") or {}).get("priority_nodes") or []
    if not nodes:
        return []
    has_q = "q_rewiring" in nodes[0]
    has_route = "selection_route" in nodes[0]
    head = "| # | Metabolite | Rewired correlations | Share of its correlations |"
    sep = "| --- | --- | --- | --- |"
    if has_q:
        head += " q (vs chance) |"
        sep += " --- |"
    if has_route:
        head += " Selected by |"
        sep += " --- |"
    head += " Abundance |"
    sep += " --- |"
    rows = [head, sep]
    for row in nodes:
        name = str(row.get("name", row.get("feature_id")))
        if name != str(row.get("feature_id")):
            name = f"{name} ({row.get('feature_id')})"
        line = (
            f"| {row.get('rank', '')} | {name} | {int(float(row.get('rewiring_score', 0)))} | "
            f"{_num(row.get('rewiring_fraction'))} |"
        )
        if has_q:
            line += f" {fmt_p(row.get('q_rewiring'))} |"
        if has_route:
            line += f" {row.get('selection_route', '')} |"
        line += f" {str(row.get('direction', '')).replace('_', ' ')} |"
        rows.append(line)
    return rows


def _rewiring_counts(section: dict[str, Any], a: str, b: str) -> str:
    s = (section.get("rewiring") or {}).get("summary") or {}

    def n(key: str) -> int:
        return int(float(s.get(key, 0) or 0))

    return (
        f"{n('n_shared_edges')} correlations in both networks, {n(f'n_{a}_only_edges')} only in {a}, "
        f"{n(f'n_{b}_only_edges')} only in {b}; {n('n_sign_switch_edges')} change sign."
    )


def export_final_report(report: dict[str, Any], outdir: Path) -> None:
    outdir = ensure_dir(Path(outdir))

    # Portable copy: relativize paths before writing user-facing reports
    report_portable = dict(report)
    report_portable["outdir"] = _relativize_path(str(report.get("outdir", "")))
    report_portable["deliverables"] = [
        _relativize_path(str(p)) for p in report.get("deliverables", [])
    ]
    save_json(report_portable, outdir / layout.REPORT_JSON)

    sections = _sections(report)
    lines: list[str] = ["# PhenoRewire results", ""]

    facts = [f"PhenoRewire {report.get('version', '')}".strip()]
    if report.get("config_path"):
        facts.append(f"config `{report['config_path']}`")
    facts.append(f"{report.get('n_samples_used', 'NA')} samples")
    n_in = report.get("n_features_input")
    n_kept = report.get("n_features_after_preprocessing", "NA")
    facts.append(f"{n_kept} of {n_in} features after preprocessing" if n_in else f"{n_kept} features after preprocessing")
    lines += [" · ".join(facts), ""]

    # ── Result ──
    lines += ["## Result", ""]
    if sections:
        for name, section, a, b in sections:
            lines.append(f"- **{_title(name, section, a, b)}:** {comparison_verdict(section, a, b)}")
    else:
        lines.append("- No analysis produced results in this run — see Warnings below.")
    lines.append("")

    warning_lines = _collect_warnings(report)
    if warning_lines:
        lines += ["## Warnings", ""]
        lines += [f"- {line}" for line in warning_lines]
        lines.append("")

    # ── Where to look ──
    lines += ["## Start here", ""]
    produced = [(n, s, a, b) for n, s, a, b in sections if s.get("rewiring")]
    if not produced:
        lines += [
            "This run produced no rewiring network. Read the Warnings above, adjust the config, "
            "and run again. The input check and the full log are in `run_info/`.",
            "",
        ]
    for name, section, a, b in produced:
        folder = section.get("folder", "")
        folder = Path(folder).name if folder else ""
        has_null = bool((section.get("rewiring") or {}).get("null"))
        lines += [f"### {_title(name, section, a, b)} — folder `{folder}/`", ""]
        lines += ["| Open | What it tells you |", "| --- | --- |"]
        if has_null:
            lines.append(
                "| `rewiring_vs_chance.csv` | Whether the two networks differ more than an arbitrary split "
                "of the same samples would. **Read it first:** if `p_empirical` > 0.05, do not read the "
                "rest as a difference between the states. |"
            )
        lines.append(
            "| `rewired_metabolites.csv` | One row per metabolite, strongest evidence first: how many of its "
            "correlations changed, whether that is more than chance (`q_rewiring`), how it was selected "
            "and which way its abundance moved. |"
        )
        lines.append(
            "| `selected_metabolites.csv` | The metabolites that entered the networks, and why "
            "(`selection_route`). |"
        )
        lines.append(
            "| `rewiring_network.graphml` | Both networks in one file, for Cytoscape: colour edges by "
            "`state`, size nodes by `rewiring_score`. |"
        )
        lines += [
            "",
            "Supporting tables are in `details/`: the statistics of every feature tested (`selection/`), "
            "each network (`network_*/`), and every correlation with its state (`rewiring/edges.csv`).",
            "",
        ]
        top = _top_rewired_rows(section)
        if top:
            lines += [f"**Most rewired metabolites** (first {len(top) - 2}, by q then rewiring score; "
                      f"all in `rewired_metabolites.csv`)", ""]
            lines += top
            lines.append("")

    lines += [
        "`run_info/` holds the run log, the input check, the configuration as it was applied, "
        "and the preprocessed data.",
        "",
    ]

    # ── Details ──
    for name, section, a, b in sections:
        lines += [f"## Details: {_title(name, section, a, b)}", ""]
        what = "abundance" if name == "Phenotype" else "association with time"
        lines.append(f"- **Selection:** {_selection_sentence(section, what)}")
        op = section.get("operating_range") or {}
        if op:
            n_a, n_b = op.get("n_per_state", ["NA", "NA"])
            lines.append(
                f"- **Sample size:** {n_a} and {n_b} samples; {op.get('status', 'NA')} the operating range "
                f"(≥ {op.get('recommended_min_n', 14)} per state recommended)."
            )
        if section.get("rewiring"):
            lines.append(f"- **Rewiring:** {_rewiring_counts(section, a, b)}")
        lines.append("")
        keys = ("network_ref", "network_case") if name == "Phenotype" else ("network_t1", "network_t2")
        if any(section.get(k) for k in keys):
            lines += _network_rows(section, keys)
            lines.append("")

    # ── Methods ──
    t_w = float(report.get("triage_topology_weight", 0.7))
    s_w = float(report.get("triage_selection_weight", 0.3))
    norm_label = NORM_LABELS.get(str(report.get("norm_method", "median")), "median normalisation")
    log_note = f"{norm_label} and log2(x + 1)" if report.get("log_transform", True) else f"{norm_label} (no log transform)"
    lines += [
        "## How the numbers are computed",
        "",
        f"- **Preprocessing:** presence filter, {log_note}.",
        "- **Selection:** a metabolite enters the networks if its abundance differs between the states "
        "(Spearman correlation with the state or with time, permutation p-values, Benjamini–Hochberg), "
        "or — with at least 14 samples per state — if its co-variation partners change.",
        "- **Networks:** Spearman correlations between the selected metabolites, kept above the |r| and "
        "q thresholds in the table above.",
        "- **rewiring_score** = correlations present in one network only + correlations that change sign. "
        "**rewiring_fraction** = the share of a metabolite's correlations present in one network only.",
        "- **Against chance:** the state labels are permuted across the pooled samples (each metabolite "
        "centred within its state first), both networks rebuilt with the same thresholds, and the "
        "state-specific correlations counted. `p_empirical` = (#null ≥ observed + 1) / (permutations + 1); "
        "`q_rewiring` is the Benjamini–Hochberg q-value per metabolite.",
        f"- **triage_score** (`details/network_*/nodes.csv`) = {t_w} × topology rank + {s_w} × (1 − q of selection).",
        "",
        "PhenoRewire generates hypotheses: a rewired metabolite changed its correlation partners; "
        "this does not show that it causes the difference between the states.",
        "",
    ]

    # ── Reproducibility ──
    lines += ["## Reproducibility", "", "```"]
    lines.append(f"version:        {report.get('version', 'NA')}")
    lines.append(f"random_seed:    {report.get('random_seed', 'NA')}")
    lines.append(f"analysis_mode:  {report.get('analysis_mode', 'NA')}")
    lines.append(f"normalisation:  {report.get('norm_method', 'median')}, log2(x + 1): {report.get('log_transform', True)}")
    for section_name, section in [("phenotype", report.get("phenotype", {})), ("temporal", report.get("temporal", {}))]:
        sel = section.get("selection", {})
        if sel:
            lines.append(
                f"{section_name} selection: FDR {sel.get('fdr_alpha', 'NA')} "
                f"(applied {sel.get('applied_fdr_alpha', 'NA')}), {sel.get('n_permutations', 'NA')} permutations"
            )
        for net_key in ["network_ref", "network_case", "network_t1", "network_t2"]:
            net = section.get(net_key, {})
            if net:
                lines.append(
                    f"network {net.get('label', net_key)}: |r| >= {net.get('applied_abs_r_thr', 'NA')}, "
                    f"q <= {net.get('applied_fdr_alpha', 'NA')}"
                )
    lines += [
        "```",
        "",
        "The same configuration (`run_info/config_used.yaml`) on the same data gives the same results.",
        "",
    ]

    (outdir / layout.REPORT_MD).write_text("\n".join(lines), encoding="utf-8")
