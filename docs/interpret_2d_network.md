# 2D Molecular Network Annotation

This document explains the optional 2D molecular network annotation workflows. `README.md` is the user-facing source of truth; this page adds more detail about the outputs.

---

## 1. Recommended workflow: package-integrated annotation tables

Use this workflow when you want practical CSV outputs for annotation review.

Run:

```bash
phenorewire --config config_phenotype.yaml --run-2dnetwork
```

This package-integrated workflow:

- takes the top PhenoRewire-prioritized rewired features
- looks for the same node IDs in an external molecular/chemical similarity network
- searches annotated neighbours within a `k`-hop neighbourhood
- traverses only edges above a similarity threshold
- writes CSV tables so you do not have to inspect the network manually

This is best described as neighbourhood-based annotation propagation over molecular networks. It is not a random walk, diffusion, or PageRank method.

---

## 2. Required inputs

The integrated workflow depends on the following inputs.

### PhenoRewire outputs

- `priority_rewired_nodes.csv`
- a rewiring GraphML produced by the same run

If `PHENOREWIRE_NETWORK` is not set in the config, the CLI automatically looks for:

- `triage/rewiring_network_phenotype.graphml`
- or `triage/rewiring_network_temporal.graphml`

### External molecular network GraphML

Set this with:

- `SPECREBOOT_NETWORK`

Requirements:

- node IDs must match the PhenoRewire feature IDs you want to expand
- edges should contain a numeric similarity attribute such as `cosine`

The current code checks edge attributes in this order:

- `cosine`
- `weight`
- `similarity`
- `score`

### Annotation node attributes

The current code looks for one of these node attributes when reporting neighbour annotations:

- `consensus_annotation`
- `annotation`
- `name`
- `label`
- `compound`

For pathway or class context, it looks for:

- `NPC#pathway`
- `npc_pathway`
- `pathway`
- `chemical_class`
- `superclass`

---

## 3. Relevant config keys

These are read from the main PhenoRewire YAML config.

| Key | Meaning |
|-----|---------|
| `SPECREBOOT_NETWORK` | Path to the molecular/chemical similarity GraphML |
| `PHENOREWIRE_NETWORK` | Optional explicit rewiring GraphML path |
| `MN_HOP_DEPTH` | Maximum neighbourhood depth |
| `MN_MIN_COSINE` | Minimum similarity required to traverse an edge |
| `ANNOTATION_TOP_N` | Number of top features to expand |

---

## 4. Main CSV outputs

The integrated workflow writes files under:

- `results/<run>/2dnetwork/`

### `annotation_summary.csv`

One row per expanded priority feature.

Important columns:

| Column | Meaning |
|--------|---------|
| `feature_id` | PhenoRewire feature ID |
| `in_phenorewire` | Whether the feature is present in the rewiring GraphML |
| `in_specreboot` | Whether the feature is present in the molecular network |
| `n_mn_neighbors` | Number of reachable molecular-network neighbours after filtering |
| `n_annotated_neighbors` | Reachable neighbours with usable annotations |
| `n_unannotated_neighbors` | Reachable neighbours without usable annotations |
| `annotation_confidence` | `n_annotated_neighbors / n_mn_neighbors` |
| `best_neighbor_cosine` | Best retained path cosine among neighbours |
| `top_annotation` | Annotation of the highest-ranked annotated neighbour |
| `rewiring_score` | Rewiring score carried through from PhenoRewire |
| `triage_score` | Triage score if available in the priority input |

### `annotation_expansion.csv`

One row per `(feature, neighbour)` pair.

Important columns:

| Column | Meaning |
|--------|---------|
| `feature_id` | Priority feature |
| `neighbor_id` | Reachable molecular-network neighbour |
| `hop_distance` | Hop distance from the priority feature |
| `edge_cosine` | Minimum edge similarity retained along the traversed path |
| `neighbor_annotation` | Neighbour annotation if present |
| `neighbor_is_annotated` | Boolean annotation flag |
| `neighbor_pathway` | Neighbour pathway/class |
| `annotation_confidence` | Feature-level confidence score copied onto the row |

### `annotation_first_neighbours.csv`

This is the fastest table to inspect manually.

It contains first-hop neighbours only and is designed for spreadsheet-style review.

Important columns:

| Column | Meaning |
|--------|---------|
| `feature_id` | Priority feature |
| `rewiring_score` | Rewiring score from PhenoRewire |
| `neighbour_id` | Direct molecular-network neighbour |
| `cosine_similarity` | Edge similarity |
| `neighbour_annotation` | Suggested neighbour annotation |
| `neighbour_pathway` | Suggested neighbour pathway/class |

---

## 5. How the current implementation works

The package-integrated workflow uses an undirected `networkx` graph and a breadth-first search.

For each priority feature:

1. start from the matching molecular-network node
2. expand neighbours up to `MN_HOP_DEPTH`
3. discard edges with similarity `< MN_MIN_COSINE`
4. collect reachable neighbours
5. count how many of those neighbours have usable annotations

The code then computes:

```text
annotation_confidence = annotated_neighbours / total_reachable_neighbours
```

So the current output is a practical neighbourhood-based annotation score, not a diffusion score.

---

## 6. Interpreting the tables

### High-confidence candidates

Look for features with:

- high `rewiring_score`
- non-zero `n_mn_neighbors`
- high `annotation_confidence`
- coherent neighbour pathway/class labels

These are strong candidates for follow-up because they are both rewired and chemically contextualized.

### Low-confidence candidates

Common causes:

- the feature is missing from the molecular network
- reachable neighbours are mostly unannotated
- the cosine threshold is too strict for the available network

### No neighbours found

Check:

- node-ID matching between the rewiring GraphML and the molecular network
- presence of a usable similarity attribute
- whether `MN_MIN_COSINE` is too high

---

## 7. Optional standalone figure workflow

The repository also contains a standalone figure-generation script:

```bash
python 2d_network/multilayer_2d_notebook.py --config 2d_network/config_phenotype.yaml
```

Use this only when you want the publication-style two-layer figure.

Install plotting extras first:

```bash
pip install "phenorewire[viz]"
```

This workflow is optional and separate from the package-integrated CSV workflow above.

---

## 8. Practical interpretation notes

- The integrated workflow does not modify the PhenoRewire rewiring scores.
- Absence from the molecular network does not mean the feature is biologically unimportant.
- A high `annotation_confidence` means the local neighbourhood is well annotated, not that the suggested annotation is proven.
