from __future__ import annotations

import argparse
import json
import logging
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib
import networkx as nx
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

logger = logging.getLogger("multilayer_2d_notebook")


DEFAULT_CONFIG: dict[str, Any] = {
    "input": {},
    "node_keys": {
        "shared_id": ["feature_id", "shared name", "id", "name", "__nodeid__"],
        "label": ["name", "label", "shared name", "feature_id", "id", "__nodeid__"],
        "chemical_label": ["consensus_annotation", "name", "label", "shared name", "feature_id", "id", "__nodeid__"],
    },
    "style_keys": {
        "edge_weight": ["weight", "correlation", "graphics_width", "width", "r_PathogenOnly", "r_3h", "r_6h"],
    },
    "focus": {
        "analysis_type": "auto",
        "top_node_limit_before_prioritization": 50,
        "keep_top_fraction": 0.5,
        "keep_top_edge_fraction": None,
        "remove_singletons": True,
        "top_layout_k": 1.25,
        "top_layout_iterations": 300,
        "top_y_offset": 3.9,
        "time_direction_attr": "time_direction",
        "score_attr": "rewiring_score",
        "state_attr": "state",
        "first_group_patterns": ["3h", "first", "time1", "t1", "PathogenOnly"],
        "second_group_patterns": ["6h", "second", "time2", "t2", "UrobiomeOnly"],
        "first_group_label": "3h",
        "second_group_label": "6h",
        "first_group_color": "#4F81BD",
        "second_group_color": "#E67E22",
        "shared_edge_color": "#8C8C8C",
        "edge_alpha": 0.40,
        "shared_edge_alpha": 0.48,
        "node_size_min": 200.0,
        "node_size_max": 600.0,
        "edge_width_min": 0.9,
        "edge_width_max": 3.2,
    },
    "comparison": {
        "blue_group": None,
        "red_group": None,
        "direction_attr": "direction",
    },
    "label_overrides": {
        "top": {},
        "bottom": {},
    },
    "chemical": {
        "annotation_id_keys": ["feature_id", "row ID", "rowID", "id"],
        "pathway_attr": "NPC#pathway",
        "annotation_label_attr": "consensus_annotation",
        "max_nodes_before_neighbor_pruning": 140,
        "matched_size": 285.0,
        "neighbor_size": 125.0,
        "neighbor_annotated_size": 160.0,
        "label_max_count": 16,
        "edge_color": "#94A3B8",
        "edge_alpha": 0.42,
        "edge_width_min": 0.6,
        "edge_width_max": 2.0,
        "palette": [
            "#F4B6C2",
            "#F6C89F",
            "#F9E79F",
            "#C7E9B4",
            "#A8DADC",
            "#BDE0FE",
            "#CDB4DB",
            "#D8E2DC",
            "#FFD6A5",
            "#E2F0CB",
            "#B8E0D2",
            "#CDE7F0",
        ],
    },
    "layout": {
        "seed": 42,
        "x_span": 10.5,
        "top_y_span": 3.7,
        "bottom_anchor_radius": 0.62,
        "bottom_anchor_radius_step": 0.18,
        "floating_radius": 5.6,
        "top_open_scale_x": 1.0,
        "top_open_scale_y": 1.22,
        "bottom_open_scale_x": 1.0,
        "bottom_open_scale_y": 1.24,
        "min_interlayer_spacing_mm": 2.0,
    },
    "labels": {
        "top_font_size": 7,
        "bottom_matched_font_size": 7,
        "bottom_neighbor_font_size": 6,
        "bottom_label_offset": 0.28,
        "top_label_offset": 0.18,
    },
    "output": {
        "basename": "multilayer_2d_notebook",
        "figure_width": 18,
        "figure_height": 13,
        "dpi": 300,
        "write_png": True,
        "write_svg": True,
        "write_pdf": True,
        "write_summary": True,
        "write_caption": True,
    },
    "titles": {
        "top": "Layer 2: phenotype network",
        "bottom": "Layer 1: chemical neighborhood",
        "show": False,
    },
    "panels": {
        "top_color": "#FCE6C9",
        "bottom_color": "#DCEAF7",
        "top_alpha": 0.14,
        "bottom_alpha": 0.18,
        "top_pad": 0.40,
        "bottom_pad": 0.50,
    },
}


@dataclass
class NodeData:
    node_id: str
    original_id: str
    label: str
    attrs: dict[str, Any]


@dataclass
class EdgeData:
    source: str
    target: str
    weight: float
    attrs: dict[str, Any]


def _deep_merge(base: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def _safe_str(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def _resolve_attr(attrs: dict[str, Any], keys: Iterable[str], node_id: str) -> Any:
    for key in keys:
        if key == "__nodeid__":
            return node_id
        if key in attrs and _safe_str(attrs[key]) != "":
            return attrs[key]
    return None


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    idx = int(round((len(ordered) - 1) * max(0.0, min(1.0, fraction))))
    return ordered[idx]


def _to_simple_graph(graph: nx.Graph) -> nx.Graph:
    simple = nx.Graph()
    for node_id, attrs in graph.nodes(data=True):
        simple.add_node(str(node_id), **dict(attrs))
    if graph.is_multigraph():
        for source, target, _key, attrs in graph.edges(data=True, keys=True):
            source_id = str(source)
            target_id = str(target)
            if source_id != target_id and not simple.has_edge(source_id, target_id):
                simple.add_edge(source_id, target_id, **dict(attrs))
    else:
        for source, target, attrs in graph.edges(data=True):
            source_id = str(source)
            target_id = str(target)
            if source_id != target_id and not simple.has_edge(source_id, target_id):
                simple.add_edge(source_id, target_id, **dict(attrs))
    return simple


def _load_graph(path: Path) -> nx.Graph:
    return _to_simple_graph(nx.read_graphml(path))


def _extract_nodes(graph: nx.Graph, key_cfg: dict[str, Any], label_key: str) -> tuple[dict[str, NodeData], dict[str, str]]:
    nodes: dict[str, NodeData] = {}
    lookup: dict[str, str] = {}
    for node_id, attrs in graph.nodes(data=True):
        node_id_str = str(node_id)
        original_id = _safe_str(_resolve_attr(attrs, key_cfg["shared_id"], node_id_str)) or node_id_str
        label = _safe_str(_resolve_attr(attrs, key_cfg[label_key], node_id_str)) or original_id
        nodes[node_id_str] = NodeData(node_id=node_id_str, original_id=original_id, label=label, attrs=dict(attrs))
        lookup[original_id] = node_id_str
    return nodes, lookup


def _apply_label_overrides(
    nodes: dict[str, NodeData],
    overrides: dict[str, str] | None,
) -> dict[str, NodeData]:
    if not overrides:
        return nodes
    updated: dict[str, NodeData] = {}
    for node_id, node in nodes.items():
        new_label = overrides.get(node.original_id, node.label)
        updated[node_id] = NodeData(
            node_id=node.node_id,
            original_id=node.original_id,
            label=new_label,
            attrs=node.attrs,
        )
    return updated


def _extract_edges(graph: nx.Graph, edge_weight_keys: list[str]) -> list[EdgeData]:
    edges: list[EdgeData] = []
    for source, target, attrs in graph.edges(data=True):
        weight = None
        for key in edge_weight_keys:
            weight = _safe_float(attrs.get(key))
            if weight is not None:
                break
        edges.append(EdgeData(str(source), str(target), abs(weight if weight is not None else 1.0), dict(attrs)))
    return edges


def _remove_singletons(graph: nx.Graph) -> nx.Graph:
    keep = [node_id for node_id in graph.nodes() if graph.degree[node_id] > 0]
    return graph.subgraph(keep).copy()


def _focus_analysis_type(nodes: dict[str, NodeData], config: dict[str, Any]) -> str:
    requested = str(config["focus"].get("analysis_type", "auto")).lower()
    if requested in {"temporal", "comparison"}:
        return requested
    if any("time_direction" in node.attrs for node in nodes.values()):
        return "temporal"
    return "comparison"


def _value_matches(value: str, patterns: Iterable[str]) -> bool:
    lower = value.lower()
    return any(pattern.lower() in lower for pattern in patterns if pattern)


def _prioritize_top_graph(
    graph: nx.Graph,
    nodes: dict[str, NodeData],
    edges: list[EdgeData],
    config: dict[str, Any],
) -> tuple[nx.Graph, dict[str, NodeData], list[EdgeData], dict[str, Any]]:
    focus_cfg = config["focus"]
    score_attr = str(focus_cfg["score_attr"])
    summary = {"prioritized": False, "removed_singletons": 0}

    working = graph.copy()
    if focus_cfg.get("remove_singletons", True):
        before = working.number_of_nodes()
        working = _remove_singletons(working)
        summary["removed_singletons"] = before - working.number_of_nodes()

    if working.number_of_nodes() <= int(focus_cfg["top_node_limit_before_prioritization"]):
        kept_nodes = {node_id: nodes[node_id] for node_id in working.nodes()}
        kept_edges = [edge for edge in edges if edge.source in kept_nodes and edge.target in kept_nodes]
        return working, kept_nodes, kept_edges, summary

    score_threshold = _percentile(
        [
            score
            for score in (_safe_float(nodes[node_id].attrs.get(score_attr)) for node_id in working.nodes())
            if score is not None
        ],
        1.0 - float(focus_cfg["keep_top_fraction"]),
    )
    if score_threshold is None:
        kept_nodes = {node_id: nodes[node_id] for node_id in working.nodes()}
        kept_edges = [edge for edge in edges if edge.source in kept_nodes and edge.target in kept_nodes]
        return working, kept_nodes, kept_edges, summary

    keep_node_ids = [
        node_id
        for node_id in working.nodes()
        if (_safe_float(nodes[node_id].attrs.get(score_attr)) or float("-inf")) >= score_threshold
    ]
    prioritized = working.subgraph(keep_node_ids).copy()
    prioritized = _remove_singletons(prioritized)

    prioritized_edges = [edge for edge in edges if edge.source in prioritized and edge.target in prioritized]
    if prioritized_edges:
        edge_fraction = focus_cfg.get("keep_top_edge_fraction")
        if edge_fraction in (None, "", 0):
            edge_fraction = float(focus_cfg["keep_top_fraction"])
        edge_threshold = _percentile([edge.weight for edge in prioritized_edges], 1.0 - float(edge_fraction))
        if edge_threshold is not None:
            prioritized_edges = [edge for edge in prioritized_edges if edge.weight >= edge_threshold]
            edge_graph = nx.Graph()
            edge_graph.add_nodes_from(prioritized.nodes(data=True))
            edge_graph.add_edges_from((edge.source, edge.target, edge.attrs) for edge in prioritized_edges)
            prioritized = _remove_singletons(edge_graph)
            prioritized_edges = [edge for edge in prioritized_edges if edge.source in prioritized and edge.target in prioritized]

    summary["prioritized"] = True
    kept_nodes = {node_id: nodes[node_id] for node_id in prioritized.nodes()}
    return prioritized, kept_nodes, prioritized_edges, summary


def _load_annotation_table(path: Path, config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    table = pd.read_csv(path, low_memory=False)
    id_col = None
    for candidate in config["chemical"]["annotation_id_keys"]:
        if candidate in table.columns:
            id_col = candidate
            break
    if id_col is None:
        raise ValueError("No feature-id column found in annotation table.")
    table[id_col] = table[id_col].astype(str).str.strip()
    table = table.drop_duplicates(subset=[id_col], keep="first")
    return table.set_index(id_col).to_dict(orient="index")


def _annotation_present(annotation: dict[str, Any], label_attr: str) -> bool:
    return _safe_str(annotation.get(label_attr, "")) != ""


def _is_feature_id_label(label: str, original_id: str) -> bool:
    clean_label = _safe_str(label)
    clean_id = _safe_str(original_id)
    if not clean_label:
        return False
    if clean_label == clean_id:
        return True
    return clean_label.isdigit()


def _build_bottom_graph(
    chemical_graph: nx.Graph,
    chemical_nodes: dict[str, NodeData],
    chemical_lookup: dict[str, str],
    matched_original_ids: set[str],
    annotation_map: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> tuple[nx.Graph, set[str], set[str]]:
    matched_node_ids = {chemical_lookup[original_id] for original_id in matched_original_ids if original_id in chemical_lookup}
    keep_node_ids = set(matched_node_ids)
    direct_neighbor_ids: set[str] = set()
    for node_id in matched_node_ids:
        direct_neighbor_ids.update(str(neighbor) for neighbor in chemical_graph.neighbors(node_id))
    direct_neighbor_ids.difference_update(matched_node_ids)

    max_nodes = int(config["chemical"]["max_nodes_before_neighbor_pruning"])
    total_nodes = len(keep_node_ids) + len(direct_neighbor_ids)
    if total_nodes > max_nodes:
        label_attr = str(config["chemical"]["annotation_label_attr"])
        ranked_neighbors = []
        for node_id in direct_neighbor_ids:
            original_id = chemical_nodes[node_id].original_id
            annotation = annotation_map.get(original_id, {})
            annotated = _annotation_present(annotation, label_attr)
            max_weight = 0.0
            for neighbor in chemical_graph.neighbors(node_id):
                edge_attrs = chemical_graph.edges[node_id, str(neighbor)]
                max_weight = max(max_weight, _safe_float(edge_attrs.get("weight")) or 0.0)
            ranked_neighbors.append((0 if annotated else 1, -max_weight, chemical_nodes[node_id].label.lower(), node_id))
        ranked_neighbors.sort()
        quota = max(max_nodes - len(keep_node_ids), 0)
        direct_neighbor_ids = {node_id for *_rest, node_id in ranked_neighbors[:quota]}

    keep_node_ids.update(direct_neighbor_ids)
    return chemical_graph.subgraph(sorted(keep_node_ids)).copy(), matched_node_ids, direct_neighbor_ids


def _normalize_positions(pos: dict[str, tuple[float, float]], x_span: float, y_span: float) -> dict[str, tuple[float, float]]:
    xs = [coord[0] for coord in pos.values()]
    ys = [coord[1] for coord in pos.values()]
    x_mid = (min(xs) + max(xs)) / 2.0
    y_mid = (min(ys) + max(ys)) / 2.0
    x_range = max(max(xs) - min(xs), 1e-9)
    y_range = max(max(ys) - min(ys), 1e-9)
    return {
        node_id: ((coord[0] - x_mid) * (x_span / x_range), (coord[1] - y_mid) * (y_span / y_range))
        for node_id, coord in pos.items()
    }


def _stretch_positions(
    pos: dict[str, tuple[float, float]],
    x_scale: float,
    y_scale: float,
) -> dict[str, tuple[float, float]]:
    return {node_id: (x * x_scale, y * y_scale) for node_id, (x, y) in pos.items()}


def _enforce_min_interlayer_spacing(
    top_pos: dict[str, tuple[float, float]],
    bottom_pos: dict[str, tuple[float, float]],
    matched_original_ids: set[str],
    top_lookup: dict[str, str],
    bottom_lookup: dict[str, str],
    bottom_anchor_map: dict[str, str],
    config: dict[str, Any],
) -> tuple[dict[str, tuple[float, float]], dict[str, tuple[float, float]]]:
    layout_cfg = config["layout"]
    min_mm = float(layout_cfg.get("min_interlayer_spacing_mm", 0.0))
    if min_mm <= 0:
        return top_pos, bottom_pos

    matched_pairs = []
    for original_id in sorted(matched_original_ids):
        top_id = top_lookup.get(original_id)
        bottom_id = bottom_lookup.get(original_id)
        if top_id in top_pos and bottom_id in bottom_pos:
            matched_pairs.append((original_id, top_id, bottom_id, float(top_pos[top_id][0])))
    if len(matched_pairs) < 2:
        return top_pos, bottom_pos

    all_x = [coord[0] for coord in top_pos.values()] + [coord[0] for coord in bottom_pos.values()]
    x_range = max(max(all_x) - min(all_x), 1e-9)
    fig_width_in = float(config["output"].get("figure_width", 14))
    plot_fraction = 0.68
    desired_gap = x_range * ((min_mm / 25.4) / max(fig_width_in * plot_fraction, 1e-9))
    if desired_gap <= 0:
        return top_pos, bottom_pos

    matched_pairs.sort(key=lambda item: item[3])
    target_x: dict[str, float] = {}
    previous_x: float | None = None
    for original_id, top_id, _bottom_id, current_x in matched_pairs:
        new_x = current_x if previous_x is None else max(current_x, previous_x + desired_gap)
        target_x[original_id] = new_x
        previous_x = new_x

    center_old = sum(item[3] for item in matched_pairs) / len(matched_pairs)
    center_new = sum(target_x[item[0]] for item in matched_pairs) / len(matched_pairs)
    center_shift = center_new - center_old
    for original_id in list(target_x):
        target_x[original_id] = target_x[original_id] - center_shift

    top_shifted = dict(top_pos)
    bottom_shifted = dict(bottom_pos)
    delta_by_anchor: dict[str, float] = {}
    for original_id, top_id, bottom_id, current_x in matched_pairs:
        delta = target_x[original_id] - current_x
        delta_by_anchor[bottom_id] = delta
        tx, ty = top_shifted[top_id]
        bx, by = bottom_shifted[bottom_id]
        top_shifted[top_id] = (tx + delta, ty)
        bottom_shifted[bottom_id] = (bx + delta, by)

    for node_id, anchor_id in bottom_anchor_map.items():
        delta = delta_by_anchor.get(anchor_id, 0.0)
        if abs(delta) <= 1e-12 or node_id not in bottom_shifted:
            continue
        x, y = bottom_shifted[node_id]
        bottom_shifted[node_id] = (x + delta, y)

    return top_shifted, bottom_shifted


def _force_vertical_interlayer_edges(
    top_pos: dict[str, tuple[float, float]],
    bottom_pos: dict[str, tuple[float, float]],
    matched_original_ids: set[str],
    top_lookup: dict[str, str],
    bottom_lookup: dict[str, str],
) -> tuple[dict[str, tuple[float, float]], dict[str, tuple[float, float]]]:
    top_aligned = dict(top_pos)
    bottom_aligned = dict(bottom_pos)
    for original_id in matched_original_ids:
        top_id = top_lookup.get(original_id)
        bottom_id = bottom_lookup.get(original_id)
        if top_id not in top_aligned or bottom_id not in bottom_aligned:
            continue
        top_x, top_y = top_aligned[top_id]
        _bottom_x, bottom_y = bottom_aligned[bottom_id]
        bottom_aligned[bottom_id] = (top_x, bottom_y)
    return top_aligned, bottom_aligned


def _build_layouts(
    top_graph: nx.Graph,
    bottom_graph: nx.Graph,
    top_lookup: dict[str, str],
    bottom_lookup: dict[str, str],
    matched_original_ids: set[str],
    bottom_nodes: dict[str, NodeData],
    config: dict[str, Any],
) -> tuple[dict[str, tuple[float, float]], dict[str, tuple[float, float]]]:
    if top_graph.number_of_nodes() == 0:
        return {}, {}
    layout_cfg = config["layout"]
    top_pos_raw = nx.spring_layout(
        top_graph,
        seed=int(layout_cfg["seed"]),
        k=float(config["focus"]["top_layout_k"]),
        iterations=int(config["focus"]["top_layout_iterations"]),
    )
    top_pos = _normalize_positions(top_pos_raw, float(layout_cfg["x_span"]), float(layout_cfg["top_y_span"]))
    top_pos = _stretch_positions(
        top_pos,
        float(layout_cfg.get("top_open_scale_x", 1.0)),
        float(layout_cfg.get("top_open_scale_y", 1.0)),
    )
    dy = float(config["focus"]["top_y_offset"])
    top_pos = {node_id: (x, y + dy) for node_id, (x, y) in top_pos.items()}

    bottom_pos: dict[str, tuple[float, float]] = {}
    bottom_anchor_map: dict[str, str] = {}
    for original_id in sorted(matched_original_ids):
        if original_id in top_lookup and original_id in bottom_lookup:
            x, y = top_pos[top_lookup[original_id]]
            bottom_pos[bottom_lookup[original_id]] = (x, y - dy)

    anchor_groups: dict[str, list[str]] = {}
    floating: list[str] = []
    for node_id in bottom_graph.nodes():
        if node_id in bottom_pos:
            continue
        matched_neighbors = [str(neighbor) for neighbor in bottom_graph.neighbors(node_id) if str(neighbor) in bottom_pos]
        if matched_neighbors:
            anchor_id = sorted(matched_neighbors)[0]
            anchor_groups.setdefault(anchor_id, []).append(node_id)
            bottom_anchor_map[node_id] = anchor_id
        else:
            floating.append(node_id)

    for anchor_node_id, node_ids in anchor_groups.items():
        cx, cy = bottom_pos[anchor_node_id]
        ordered = sorted(node_ids, key=lambda n: bottom_nodes[n].label.lower())
        for idx, node_id in enumerate(ordered):
            angle = 2.0 * math.pi * idx / max(len(ordered), 1)
            radius = float(layout_cfg["bottom_anchor_radius"]) + float(layout_cfg["bottom_anchor_radius_step"]) * (idx % 4)
            bottom_pos[node_id] = (cx + radius * math.cos(angle), cy + radius * math.sin(angle))

    if floating:
        radius = float(layout_cfg["floating_radius"])
        for idx, node_id in enumerate(sorted(floating, key=lambda n: bottom_nodes[n].label.lower())):
            angle = 2.0 * math.pi * idx / max(len(floating), 1)
            bottom_pos[node_id] = (radius * math.cos(angle), radius * math.sin(angle) - 1.0)

    bottom_pos = _stretch_positions(
        bottom_pos,
        float(layout_cfg.get("bottom_open_scale_x", 1.0)),
        float(layout_cfg.get("bottom_open_scale_y", 1.0)),
    )
    top_pos, bottom_pos = _enforce_min_interlayer_spacing(
        top_pos,
        bottom_pos,
        matched_original_ids,
        top_lookup,
        bottom_lookup,
        bottom_anchor_map,
        config,
    )
    top_pos, bottom_pos = _force_vertical_interlayer_edges(
        top_pos,
        bottom_pos,
        matched_original_ids,
        top_lookup,
        bottom_lookup,
    )
    return top_pos, bottom_pos


def _scale(values: dict[str, float], min_size: float, max_size: float) -> dict[str, float]:
    if not values:
        return {}
    vmin = min(values.values())
    vmax = max(values.values())
    if abs(vmax - vmin) < 1e-9:
        return {key: (min_size + max_size) / 2.0 for key in values}
    return {
        key: min_size + ((value - vmin) / (vmax - vmin)) * (max_size - min_size)
        for key, value in values.items()
    }


def _group_for_top_node(node: NodeData, analysis_type: str, config: dict[str, Any]) -> str:
    focus_cfg = config["focus"]
    if analysis_type == "temporal":
        value = _safe_str(node.attrs.get(focus_cfg["time_direction_attr"]))
        if _value_matches(value, focus_cfg["first_group_patterns"]):
            return "first"
        if _value_matches(value, focus_cfg["second_group_patterns"]):
            return "second"
    else:
        value = _safe_str(node.attrs.get(config["comparison"].get("direction_attr", "direction")))
        blue_group = _safe_str(config["comparison"].get("blue_group"))
        red_group = _safe_str(config["comparison"].get("red_group"))
        if blue_group and blue_group.lower() in value.lower():
            return "first"
        if red_group and red_group.lower() in value.lower():
            return "second"
    return "other"


def _top_edge_style(edge: EdgeData, analysis_type: str, config: dict[str, Any]) -> tuple[str, float]:
    focus_cfg = config["focus"]
    state = _safe_str(edge.attrs.get(focus_cfg["state_attr"]))
    if "shared" in state.lower():
        return str(focus_cfg["shared_edge_color"]), float(focus_cfg["shared_edge_alpha"])
    if analysis_type == "temporal":
        if _value_matches(state, focus_cfg["first_group_patterns"]):
            return str(focus_cfg["first_group_color"]), float(focus_cfg["edge_alpha"])
        if _value_matches(state, focus_cfg["second_group_patterns"]):
            return str(focus_cfg["second_group_color"]), float(focus_cfg["edge_alpha"])
    if _value_matches(state, [_safe_str(config["comparison"].get("blue_group"))]):
        return str(focus_cfg["first_group_color"]), float(focus_cfg["edge_alpha"])
    if _value_matches(state, [_safe_str(config["comparison"].get("red_group"))]):
        return str(focus_cfg["second_group_color"]), float(focus_cfg["edge_alpha"])
    return str(focus_cfg["shared_edge_color"]), float(focus_cfg["shared_edge_alpha"])


def _chemical_palette(bottom_nodes: dict[str, NodeData], annotation_map: dict[str, dict[str, Any]], config: dict[str, Any]) -> dict[str, str]:
    pathway_attr = str(config["chemical"]["pathway_attr"])
    values = []
    for node in bottom_nodes.values():
        values.append(_safe_str(annotation_map.get(node.original_id, {}).get(pathway_attr)) or "Unassigned")
    unique = sorted(set(values))
    palette = config["chemical"]["palette"]
    return {value: palette[idx % len(palette)] for idx, value in enumerate(unique)}


def _safe_save(fig: plt.Figure, path: Path, dpi: int) -> str:
    try:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        return str(path)
    except PermissionError:
        fallback = path.with_name(f"{path.stem}_new{path.suffix}")
        fig.savefig(fallback, dpi=dpi, bbox_inches="tight")
        return str(fallback)


def _bottom_label_xy(x: float, y: float, offset: float, cx: float = 0.0, cy: float = 0.0) -> tuple[float, float, str, str]:
    """Place label radially away from the cluster centroid (cx, cy)."""
    dx = x - cx
    dy = y - cy
    norm = math.sqrt(dx * dx + dy * dy)
    if norm < 1e-9:
        dx, dy = 1.0, 0.0
    else:
        dx /= norm
        dy /= norm
    lx = x + dx * offset * 1.5
    ly = y + dy * offset * 1.5
    ha = "left" if dx >= 0 else "right"
    va = "bottom" if dy >= 0 else "top"
    return (lx, ly, ha, va)


def build_notebook_multilayer(
    top_graphml: Path,
    chemical_graphml: Path,
    annotation_table: Path,
    outdir: Path,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    merged = _deep_merge(DEFAULT_CONFIG, config or {})
    outdir.mkdir(parents=True, exist_ok=True)

    top_graph = _load_graph(Path(top_graphml))
    chemical_graph = _load_graph(Path(chemical_graphml))
    annotation_map = _load_annotation_table(Path(annotation_table), merged)

    top_nodes, top_lookup = _extract_nodes(top_graph, merged["node_keys"], "label")
    top_nodes = _apply_label_overrides(top_nodes, merged.get("label_overrides", {}).get("top", {}))
    top_edges = _extract_edges(top_graph, merged["style_keys"]["edge_weight"])
    top_graph, top_nodes, top_edges, prioritization = _prioritize_top_graph(top_graph, top_nodes, top_edges, merged)
    top_lookup = {node.original_id: node_id for node_id, node in top_nodes.items()}

    if top_graph.number_of_nodes() == 0:
        logger.warning("No nodes remain in the top layer after prioritization - check remove_singletons and keep_top_fraction settings.")
        return {
            "top_nodes": 0, "top_edges": 0, "bottom_nodes": 0, "bottom_edges": 0,
            "matched_nodes": 0, "direct_neighbors": 0,
            "prioritized_top": prioritization["prioritized"],
            "removed_singletons": prioritization["removed_singletons"],
        }

    chemical_nodes, chemical_lookup = _extract_nodes(chemical_graph, merged["node_keys"], "chemical_label")
    chemical_edges = _extract_edges(chemical_graph, merged["style_keys"]["edge_weight"])

    matched_original_ids = set(top_lookup.keys()) & set(chemical_lookup.keys())
    bottom_graph, matched_bottom_node_ids, direct_neighbor_ids = _build_bottom_graph(
        chemical_graph,
        chemical_nodes,
        chemical_lookup,
        matched_original_ids,
        annotation_map,
        merged,
    )
    bottom_nodes = {node_id: chemical_nodes[node_id] for node_id in bottom_graph.nodes()}
    bottom_edges = [edge for edge in chemical_edges if edge.source in bottom_nodes and edge.target in bottom_nodes]

    top_pos, bottom_pos = _build_layouts(
        top_graph,
        bottom_graph,
        top_lookup,
        chemical_lookup,
        matched_original_ids,
        bottom_nodes,
        merged,
    )

    analysis_type = _focus_analysis_type(top_nodes, merged)
    top_score_attr = str(merged["focus"]["score_attr"])
    top_scores = {node_id: (_safe_float(node.attrs.get(top_score_attr)) or 0.0) for node_id, node in top_nodes.items()}
    top_node_sizes = _scale(top_scores, float(merged["focus"]["node_size_min"]), float(merged["focus"]["node_size_max"]))
    top_edge_sizes = _scale({f"{edge.source}|{edge.target}": edge.weight for edge in top_edges}, float(merged["focus"]["edge_width_min"]), float(merged["focus"]["edge_width_max"]))
    bottom_edge_sizes = _scale({f"{edge.source}|{edge.target}": edge.weight for edge in bottom_edges}, float(merged["chemical"]["edge_width_min"]), float(merged["chemical"]["edge_width_max"]))
    pathway_colors = _chemical_palette(bottom_nodes, annotation_map, merged)

    fig, ax = plt.subplots(figsize=(merged["output"]["figure_width"], merged["output"]["figure_height"]))
    ax.set_aspect("equal")
    ax.axis("off")

    if bottom_pos:
        bottom_x = [coord[0] for coord in bottom_pos.values()]
        bottom_y = [coord[1] for coord in bottom_pos.values()]
        ax.add_patch(Rectangle((min(bottom_x) - float(merged["panels"]["bottom_pad"]), min(bottom_y) - float(merged["panels"]["bottom_pad"])), (max(bottom_x) - min(bottom_x)) + 2 * float(merged["panels"]["bottom_pad"]), (max(bottom_y) - min(bottom_y)) + 2 * float(merged["panels"]["bottom_pad"]), facecolor=merged["panels"]["bottom_color"], edgecolor="none", alpha=float(merged["panels"]["bottom_alpha"]), zorder=0))
        if merged["titles"].get("show", False):
            ax.text(min(bottom_x) - 0.15, min(bottom_y) - 0.42, merged["titles"]["bottom"], fontsize=12, fontweight="bold")
    if top_pos:
        top_x_vals = [coord[0] for coord in top_pos.values()]
        top_y_vals = [coord[1] for coord in top_pos.values()]
        ax.add_patch(Rectangle((min(top_x_vals) - float(merged["panels"]["top_pad"]), min(top_y_vals) - float(merged["panels"]["top_pad"])), (max(top_x_vals) - min(top_x_vals)) + 2 * float(merged["panels"]["top_pad"]), (max(top_y_vals) - min(top_y_vals)) + 2 * float(merged["panels"]["top_pad"]), facecolor=merged["panels"]["top_color"], edgecolor="none", alpha=float(merged["panels"]["top_alpha"]), zorder=0))
        if merged["titles"].get("show", False):
            ax.text(min(top_x_vals) - 0.15, max(top_y_vals) + 0.28, merged["titles"]["top"], fontsize=12, fontweight="bold")

    for original_id in sorted(matched_original_ids):
        if original_id in top_lookup and original_id in chemical_lookup:
            top_id = top_lookup[original_id]
            chem_id = chemical_lookup[original_id]
            if top_id in top_pos and chem_id in bottom_pos:
                x1, y1 = bottom_pos[chem_id]
                x2, y2 = top_pos[top_id]
                ax.plot([x1, x2], [y1, y2], linestyle=(0, (2, 3)), color="#5F5F5F", alpha=0.52, linewidth=1.5, zorder=1)

    if bottom_edges:
        segments = [[bottom_pos[edge.source], bottom_pos[edge.target]] for edge in bottom_edges]
        widths = [bottom_edge_sizes[f"{edge.source}|{edge.target}"] for edge in bottom_edges]
        ax.add_collection(LineCollection(segments, colors=[merged["chemical"]["edge_color"]] * len(segments), linewidths=widths, alpha=float(merged["chemical"]["edge_alpha"]), zorder=2))

    if top_edges:
        top_edge_groups: dict[tuple[str, float], list[tuple[list, float]]] = defaultdict(list)
        for edge in top_edges:
            color, alpha = _top_edge_style(edge, analysis_type, merged)
            width = top_edge_sizes[f"{edge.source}|{edge.target}"]
            top_edge_groups[(color, alpha)].append(([top_pos[edge.source], top_pos[edge.target]], width))
        for (color, alpha), seg_widths in top_edge_groups.items():
            segs = [sw[0] for sw in seg_widths]
            widths = [sw[1] for sw in seg_widths]
            ax.add_collection(LineCollection(segs, colors=[color] * len(segs), linewidths=widths, alpha=alpha, zorder=3))

    pathway_attr = str(merged["chemical"]["pathway_attr"])
    for node_id, node in bottom_nodes.items():
        x, y = bottom_pos[node_id]
        annotation = annotation_map.get(node.original_id, {})
        pathway = _safe_str(annotation.get(pathway_attr)) or "Unassigned"
        color = pathway_colors[pathway]
        if node_id in matched_bottom_node_ids:
            size = float(merged["chemical"]["matched_size"])
            edgecolor = "black"
            linewidth = 1.4
        else:
            size = float(merged["chemical"]["neighbor_annotated_size"]) if _annotation_present(annotation, merged["chemical"]["annotation_label_attr"]) else float(merged["chemical"]["neighbor_size"])
            edgecolor = "white"
            linewidth = 0.8
        ax.scatter([x], [y], s=[size], c=[color], edgecolors=[edgecolor], linewidths=[linewidth], alpha=0.95, zorder=4)

    for node_id, node in top_nodes.items():
        x, y = top_pos[node_id]
        group = _group_for_top_node(node, analysis_type, merged)
        if group == "first":
            color = merged["focus"]["first_group_color"]
        elif group == "second":
            color = merged["focus"]["second_group_color"]
        else:
            color = "#BDBDBD"
        ax.scatter([x], [y], s=[top_node_sizes[node_id]], c=[color], edgecolors=["black"], linewidths=[1.1], alpha=0.98, zorder=5)
        if not _is_feature_id_label(node.label, node.original_id):
            ax.text(
                x,
                y,
                node.label,
                fontsize=int(merged["labels"]["top_font_size"]),
                ha="center",
                va="center",
                zorder=6,
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white", alpha=0.95)],
            )

    label_attr = str(merged["chemical"]["annotation_label_attr"])
    bottom_label_offset = float(merged["labels"]["bottom_label_offset"])

    # Cluster centroid for radial label placement (direction away from centre, not from origin)
    if bottom_pos:
        _bxs = [c[0] for c in bottom_pos.values()]
        _bys = [c[1] for c in bottom_pos.values()]
        bottom_cx = sum(_bxs) / len(_bxs)
        bottom_cy = sum(_bys) / len(_bys)
    else:
        bottom_cx = bottom_cy = 0.0

    # Max edge weight per bottom node (used to rank neighbor labels)
    node_max_weight: dict[str, float] = {}
    for edge in bottom_edges:
        node_max_weight[edge.source] = max(node_max_weight.get(edge.source, 0.0), edge.weight)
        node_max_weight[edge.target] = max(node_max_weight.get(edge.target, 0.0), edge.weight)

    # Labels for matched bottom nodes — always shown (these are the annotation-propagation targets)
    for node_id in sorted(matched_bottom_node_ids):
        node = bottom_nodes.get(node_id)
        if node is None or node_id not in bottom_pos:
            continue
        annotation = annotation_map.get(node.original_id, {})
        label = _safe_str(annotation.get(label_attr))
        if not label or _is_feature_id_label(label, node.original_id):
            continue
        x, y = bottom_pos[node_id]
        label_x, label_y, ha, va = _bottom_label_xy(x, y, bottom_label_offset, bottom_cx, bottom_cy)
        ax.text(
            label_x, label_y, label,
            fontsize=int(merged["labels"]["bottom_matched_font_size"]),
            ha=ha, va=va, zorder=6,
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white", alpha=0.95)],
        )

    # Labels for direct neighbor nodes — sorted by max edge weight, best-connected first
    neighbor_label_candidates: list[tuple[float, float, str, int, float]] = []
    for node_id, node in bottom_nodes.items():
        if node_id not in direct_neighbor_ids or node_id not in bottom_pos:
            continue
        annotation = annotation_map.get(node.original_id, {})
        label = _safe_str(annotation.get(label_attr))
        if not label or _is_feature_id_label(label, node.original_id):
            continue
        x, y = bottom_pos[node_id]
        max_w = node_max_weight.get(node_id, 0.0)
        neighbor_label_candidates.append((x, y, label, int(merged["labels"]["bottom_neighbor_font_size"]), max_w))
    neighbor_label_candidates.sort(key=lambda item: -item[4])
    for x, y, label, font_size, _ in neighbor_label_candidates[: int(merged["chemical"]["label_max_count"])]:
        label_x, label_y, ha, va = _bottom_label_xy(x, y, bottom_label_offset, bottom_cx, bottom_cy)
        ax.text(
            label_x, label_y, label,
            fontsize=font_size, ha=ha, va=va, zorder=6,
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white", alpha=0.95)],
        )

    phenotype_handles = [
        Patch(facecolor=merged["focus"]["first_group_color"], edgecolor="black", label=str(merged["focus"]["first_group_label"])),
        Patch(facecolor=merged["focus"]["second_group_color"], edgecolor="black", label=str(merged["focus"]["second_group_label"])),
        Line2D([0], [0], color=str(merged["focus"]["shared_edge_color"]), linewidth=1.8,
               alpha=float(merged["focus"]["shared_edge_alpha"]), label="Shared edges"),
    ]
    fig.subplots_adjust(right=0.74)
    legend_ax = fig.add_axes([0.77, 0.16, 0.18, 0.76])
    legend_ax.axis("off")
    phenotype_legend = legend_ax.legend(handles=phenotype_handles, title="Phenotype", loc="upper left", frameon=False, fontsize=9, title_fontsize=10)
    legend_ax.add_artist(phenotype_legend)

    pathway_handles = [Patch(facecolor=color, edgecolor="#555555", label=pathway) for pathway, color in sorted(pathway_colors.items())]
    chemical_legend = legend_ax.legend(handles=pathway_handles, title="NPC pathway", loc="lower left", frameon=False, fontsize=9, title_fontsize=10)
    legend_ax.add_artist(chemical_legend)

    basename = str(merged["output"]["basename"])
    outputs: dict[str, str] = {}
    for ext, enabled in [("png", merged["output"]["write_png"]), ("svg", merged["output"]["write_svg"]), ("pdf", merged["output"]["write_pdf"])]:
        if enabled:
            outputs[ext] = _safe_save(fig, outdir / f"{basename}.{ext}", int(merged["output"]["dpi"]))
    plt.close(fig)

    summary = {
        "top_nodes": top_graph.number_of_nodes(),
        "top_edges": top_graph.number_of_edges(),
        "bottom_nodes": bottom_graph.number_of_nodes(),
        "bottom_edges": bottom_graph.number_of_edges(),
        "matched_nodes": len(matched_original_ids),
        "direct_neighbors": len(direct_neighbor_ids),
        "prioritized_top": prioritization["prioritized"],
        "removed_singletons": prioritization["removed_singletons"],
    }
    if merged["output"].get("write_summary", True):
        summary_path = outdir / f"{basename}_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        outputs["summary_json"] = str(summary_path)
    if merged["output"].get("write_caption", True):
        first_label = str(merged["focus"]["first_group_label"])
        second_label = str(merged["focus"]["second_group_label"])
        mode_phrase = (
            f"temporal ({first_label} vs {second_label})"
            if analysis_type == "temporal"
            else f"phenotype comparison ({first_label} vs {second_label})"
        )
        caption = (
            f"Two-layer 2D network integrating a {mode_phrase} rewiring network (top layer) "
            f"with its chemical similarity neighborhood (bottom layer). "
            f"Top-layer nodes are colored by biological state "
            f"({first_label}: blue; {second_label}: orange) and scaled by rewiring score. "
            "Bottom-layer nodes are colored by NPC pathway; matched metabolites carry a black outline "
            "and their direct spectral neighbors are shown as smaller nodes. "
            "Dashed vertical connections link matched metabolites across layers. "
            "The bottom layer supports annotation propagation by exposing the chemical context "
            "of PhenoRewire-prioritized features."
        )
        caption_path = outdir / f"{basename}_caption.txt"
        caption_path.write_text(caption, encoding="utf-8")
        outputs["caption_txt"] = str(caption_path)
    summary["outputs"] = outputs
    return summary


def load_yaml_config(path: Path | str) -> dict[str, Any]:
    config_path = Path(path)
    data = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError("Config YAML must be a mapping.")
    input_cfg = data.get("input", {})
    for key in ("top_graphml", "chemical_graphml", "annotation_table", "outdir"):
        value = input_cfg.get(key)
        if value is None:
            continue
        candidate = Path(str(value))
        if not candidate.is_absolute():
            input_cfg[key] = str((config_path.parent / candidate).resolve())
    return data


_REQUIRED_INPUT_KEYS = ("top_graphml", "chemical_graphml", "annotation_table", "outdir")


def run_from_yaml(config_path: Path | str) -> dict[str, Any]:
    data = load_yaml_config(config_path)
    input_cfg = data.get("input", {})

    missing = [key for key in _REQUIRED_INPUT_KEYS if not input_cfg.get(key)]
    if missing:
        # The most common cause is handing this script a pipeline config: the two
        # config formats look similar but are not interchangeable.
        looks_like_pipeline = any(key in data for key in ("DATA_MATRIX", "GROUP_DEFINITION", "OUTDIR"))
        hint = (
            "\n\nThis looks like a PhenoRewire *pipeline* config (it has DATA_MATRIX / "
            "GROUP_DEFINITION / OUTDIR). This script takes a *figure* config instead: "
            "see 2d_network/2d_figure_phenotype.yaml for a template."
            if looks_like_pipeline
            else "\n\nSee 2d_network/2d_figure_phenotype.yaml for a template."
        )
        raise ValueError(
            f"Config '{config_path}' is missing required input key(s): {missing}.\n"
            f"Expected a top-level 'input:' block with {list(_REQUIRED_INPUT_KEYS)}."
            + hint
        )

    config = dict(data)
    config.pop("input", None)
    return build_notebook_multilayer(
        top_graphml=Path(input_cfg["top_graphml"]),
        chemical_graphml=Path(input_cfg["chemical_graphml"]),
        annotation_table=Path(input_cfg["annotation_table"]),
        outdir=Path(input_cfg["outdir"]),
        config=config,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Notebook-inspired multilayer 2D network builder.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, str(args.log_level).upper(), logging.INFO), format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    run_from_yaml(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
