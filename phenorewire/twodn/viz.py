"""
Visualization helpers for PhenoRewire network outputs.

Requires matplotlib (hard dependency of this module) and pyvis (optional, for the
interactive HTML plot).  matplotlib is imported at module level on purpose: callers
guard ``import phenorewire.twodn.viz`` with ``except ImportError`` to tell the user
to install the ``viz`` extra, and a lazy import inside the function would defeat that.
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def plot_rewiring_network(
    graphml_path: Path,
    *,
    triage_df=None,
    rewiring_df=None,
    output_path: Path,
    top_n: int = 40,
) -> None:
    """
    Produce PNG and (optionally) interactive HTML plots of the PhenoRewire rewiring network.

    Parameters
    ----------
    graphml_path : path to rewiring_network.graphml
    triage_df    : DataFrame with feature_id and triage_score (for node colouring)
    rewiring_df  : DataFrame with feature_id and rewiring_score (for node size)
    output_path  : stem path (without extension); .png and .html are appended
    top_n        : max nodes to display for readability
    """
    graphml_path = Path(graphml_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not graphml_path.exists():
        logger.warning("rewiring_network.graphml not found: %s", graphml_path)
        return

    G = nx.read_graphml(str(graphml_path))
    if G.is_directed():
        G = G.to_undirected()

    if G.number_of_nodes() == 0:
        logger.warning("Rewiring network is empty — skipping plot.")
        return

    # Optionally subset to top_n nodes by rewiring_score
    score_map: dict[str, float] = {}
    if rewiring_df is not None and not rewiring_df.empty:
        id_col = _id_col(rewiring_df)
        if id_col and "rewiring_score" in rewiring_df.columns:
            score_map = dict(zip(
                rewiring_df[id_col].astype(str),
                rewiring_df["rewiring_score"].astype(float),
            ))

    if score_map:
        sorted_nodes = sorted(G.nodes, key=lambda n: score_map.get(n, 0.0), reverse=True)
        keep = set(sorted_nodes[:top_n])
        G = G.subgraph(keep).copy()

    triage_map: dict[str, float] = {}
    if triage_df is not None and not triage_df.empty:
        id_col = _id_col(triage_df)
        if id_col and "triage_score" in triage_df.columns:
            triage_map = dict(zip(
                triage_df[id_col].astype(str),
                triage_df["triage_score"].astype(float),
            ))

    nodes = list(G.nodes)
    triage_scores = np.array([triage_map.get(n, 0.0) for n in nodes])
    rew_scores = np.array([score_map.get(n, 0.0) for n in nodes])

    # Node colours from triage_score (viridis colormap)
    cmap = plt.cm.viridis
    norm_t = mcolors.Normalize(vmin=triage_scores.min(), vmax=max(triage_scores.max(), 1e-9))
    node_colors = [cmap(norm_t(s)) for s in triage_scores]

    # Node sizes from rewiring_score (100–500)
    r_min, r_max = rew_scores.min(), rew_scores.max()
    if r_max > r_min:
        node_sizes = 100 + 400 * (rew_scores - r_min) / (r_max - r_min)
    else:
        node_sizes = np.full(len(nodes), 200.0)

    # Edge styling: colour by state attribute if present
    edge_colors = []
    for u, v, d in G.edges(data=True):
        state = str(d.get("state", "shared")).lower()
        if "only_a" in state or "ref" in state:
            edge_colors.append("#4F81BD")
        elif "only_b" in state or "case" in state:
            edge_colors.append("#E67E22")
        elif "sign_switch" in state:
            edge_colors.append("#E74C3C")
        else:
            edge_colors.append("#AAAAAA")

    pos = nx.spring_layout(G, seed=42, k=1.5)

    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors, alpha=0.5, width=1.2)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=node_sizes, alpha=0.9)
    nx.draw_networkx_labels(
        G, pos, ax=ax,
        labels={n: n for n in nodes},
        font_size=6, font_color="black",
    )
    ax.set_title("PhenoRewire Rewiring Network", fontsize=14)
    ax.axis("off")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm_t)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Triage score", shrink=0.7)

    png_path = output_path.with_suffix(".png")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved PNG: %s", png_path)

    # Interactive HTML via pyvis (optional extra)
    try:
        _plot_interactive(G, score_map, triage_map, output_path.with_suffix(".html"))
    except ImportError:
        logger.debug("pyvis not installed — skipping interactive HTML plot.")


def _plot_interactive(G, score_map: dict, triage_map: dict, html_path: Path) -> None:
    from pyvis.network import Network

    nt = Network(height="750px", width="100%", bgcolor="#1a1a2e", font_color="white")
    nt.set_options('{"physics": {"stabilization": {"iterations": 150}}}')

    for node in G.nodes:
        triage = triage_map.get(node, 0.0)
        rew = score_map.get(node, 0.0)
        size = 10 + 25 * triage
        r = int(255 * triage)
        b = int(255 * (1 - triage))
        color = f"#{r:02x}88{b:02x}"
        nt.add_node(
            node, label=node, size=size,
            color=color, title=f"triage={triage:.3f} rew={rew:.3f}",
        )

    for u, v, d in G.edges(data=True):
        state = str(d.get("state", "shared")).lower()
        if "sign_switch" in state:
            ecolor = "#E74C3C"
        elif "only" in state:
            ecolor = "#F39C12"
        else:
            ecolor = "#7F8C8D"
        nt.add_edge(u, v, color=ecolor)

    nt.save_graph(str(html_path))
    logger.info("Saved interactive HTML: %s", html_path)


def _id_col(df) -> str | None:
    for c in ("feature_id", "node_id", "id", "feature"):
        if c in df.columns:
            return c
    return df.columns[0] if len(df.columns) > 0 else None
