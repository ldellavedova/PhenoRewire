# ── PhenoRewire · Step 4 — Triage & ranking ──────────────────────────────────
# Ranks the nodes of each network by triage_score, and the metabolites of a
# comparison by how much (and how significantly) they rewire.
# triage_score = 0.7 * topology_priority_score + 0.3 * selection_support_score
#   topology_priority_score = mean percentile rank of [degree_centrality,
#     betweenness_centrality, eigenvector_centrality, weighted_degree_abs]
#   selection_support_score = 1 - q_value (configurable weights via YAML)
# rewiring_score = state_A_only_edges + state_B_only_edges + sign_switch_edges
# rewiring_fraction = (state_A_only + state_B_only) / total_union_degree  [0–1]
# Both functions return tables; run.py decides where they are written.
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def _rank_pct(series: pd.Series) -> pd.Series:
    if series.empty:
        return pd.Series(dtype=float, index=series.index)
    return series.rank(method="average", pct=True).fillna(0.0)


def rank_network_nodes(
    nodes_df: pd.DataFrame,
    *,
    topology_weight: float = 0.7,
    selection_weight: float = 0.3,
) -> pd.DataFrame:
    """
    The nodes of one network, with triage scores, best first.

    triage_score = topology_weight * topology_priority_score
                 + selection_weight * selection_support_score

    topology_priority_score = mean percentile rank of
        [degree_centrality, betweenness_centrality, eigenvector_centrality, weighted_degree_abs]
    selection_support_score = 1 - q_value (phenotype/time association strength, 0–1)

    The default weights (0.7 topology / 0.3 selection) were validated in
    dev/simulate/weight_sensitivity.py across [0.1/0.9 ... 0.9/0.1] on a synthetic dataset
    with hub-only and selection-only ground-truth features.  AUROC for ground-truth node
    recovery varies ~0.19 across the range (highest at selection-heavy 0.1/0.9, lowest at
    topology-heavy 0.9/0.1), confirming that weights matter.  The 0.7/0.3 default balances
    hub detection with phenotype-association signal.  Both weights are configurable via
    TRIAGE_TOPOLOGY_WEIGHT / TRIAGE_SELECTION_WEIGHT in the YAML config.
    """
    if nodes_df is None or nodes_df.empty:
        return pd.DataFrame(columns=["feature_id", "triage_score"])

    ranked = nodes_df.copy()
    if "feature_id" not in ranked.columns:
        ranked = ranked.reset_index(drop=False)
        ranked = ranked.rename(columns={ranked.columns[0]: "feature_id"})
    ranked = ranked.reset_index(drop=True)

    for col in ["degree_centrality", "betweenness_centrality", "eigenvector_centrality", "weighted_degree_abs"]:
        if col not in ranked.columns:
            ranked[col] = 0.0

    ranked["topology_priority_score"] = (
        _rank_pct(ranked["degree_centrality"].astype(float))
        + _rank_pct(ranked["betweenness_centrality"].astype(float))
        + _rank_pct(ranked["eigenvector_centrality"].astype(float))
        + _rank_pct(ranked["weighted_degree_abs"].astype(float))
    ) / 4.0

    for q_col in ("q_selection", "q_pheno", "q_time"):
        # q_selection is set when connectivity-based selection ran: the q-value of
        # the route that selected the feature (the smaller one when both did).
        if q_col in ranked.columns:
            q = pd.to_numeric(ranked[q_col], errors="coerce")
            ranked["selection_support_score"] = 1.0 - q.clip(lower=0.0, upper=1.0).fillna(1.0)
            break
    else:
        ranked["selection_support_score"] = 0.0

    ranked["triage_score"] = (
        float(topology_weight) * ranked["topology_priority_score"]
        + float(selection_weight) * ranked["selection_support_score"]
    )
    return ranked.sort_values(
        ["triage_score", "topology_priority_score", "feature_id"], ascending=[False, False, True]
    ).reset_index(drop=True)


def rewired_metabolites(
    selected: pd.DataFrame,
    rewiring_nodes: pd.DataFrame,
    node_null: Optional[pd.DataFrame],
    *,
    label_a: str,
    label_b: str,
    direction_col: str,
) -> pd.DataFrame:
    """
    One row per selected metabolite: who it is, how it was selected, which way it
    moved, how many of its connections changed, and whether that is more than chance.

    Metabolites without any edge in either network are kept with zero counts, so
    the table always lists every metabolite that entered the comparison.  Sorted by
    q_rewiring (when the null model ran), then by rewiring_score.
    """
    ids = selected.index.astype(str)
    out = pd.DataFrame({"feature_id": ids})
    out["name"] = selected["name"].astype(str).values if "name" in selected.columns else ids
    if "selection_route" in selected.columns:
        out["selection_route"] = selected["selection_route"].values
    if direction_col in selected.columns:
        out["direction"] = selected[direction_col].values

    count_cols = [f"{label_a}_only_edges", f"{label_b}_only_edges", "shared_edges", "sign_switch_edges"]
    rew = rewiring_nodes.copy() if rewiring_nodes is not None else pd.DataFrame(columns=["feature_id"])
    rew["feature_id"] = rew["feature_id"].astype(str)
    keep = ["feature_id", "rewiring_score", "rewiring_fraction"] + count_cols
    out = out.merge(rew[[c for c in keep if c in rew.columns]], on="feature_id", how="left")
    for col in ["rewiring_score", "rewiring_fraction"] + count_cols:
        if col not in out.columns:
            out[col] = 0
        out[col] = out[col].fillna(0)
    for col in count_cols + ["rewiring_score"]:
        out[col] = out[col].astype(int)

    sort_cols, ascending = ["rewiring_score"], [False]
    if node_null is not None and not node_null.empty:
        null = node_null.rename(
            columns={
                "null_mean_state_specific_edges": "expected_by_chance",
                "p_empirical": "p_rewiring",
                "q_bh": "q_rewiring",
            }
        )
        null["feature_id"] = null["feature_id"].astype(str)
        out = out.merge(
            null[["feature_id", "expected_by_chance", "p_rewiring", "q_rewiring"]], on="feature_id", how="left"
        )
        out["more_than_chance"] = out["q_rewiring"].astype(float) <= 0.05
        sort_cols, ascending = ["q_rewiring", "rewiring_score"], [True, False]

    out = out.sort_values(sort_cols + ["feature_id"], ascending=ascending + [True]).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    # The verdict first, then how it was selected, then the counts behind it.
    order = [
        "rank", "feature_id", "name", "rewiring_score", "rewiring_fraction", "q_rewiring",
        "more_than_chance", "selection_route", "direction",
    ] + count_cols + ["expected_by_chance", "p_rewiring"]
    return out[[c for c in order if c in out.columns]]
