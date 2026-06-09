# How to interpret PhenoRewire results — phenotype mode

This guide is written for metabolomics researchers who may have no bioinformatics background. It explains what phenotype mode does, which files to open first, and what each column means.

---

## What phenotype mode does

Phenotype mode compares two biological groups — for example, case vs reference, treated vs control, or responders vs non-responders — measured at the same timepoint. PhenoRewire first identifies which metabolite features are statistically associated with group identity. It then builds one correlation network per group using only those selected features. Finally, it compares the two networks edge by edge to quantify which metabolites changed not just in abundance but in how they co-vary with other metabolites. The output describes which features are most reorganised in their metabolic co-regulation context between the two states.

---

## Start here — your two primary output files

After every phenotype run, open these two files first:

**1. `triage/rewiring_network_phenotype.graphml`** — open in [Cytoscape](https://cytoscape.org) (free download)

Nodes are metabolite features. Node size is proportional to `rewiring_score` — larger nodes changed more. Edge color distinguishes:
- **Grey edges**: shared between both groups (conserved co-regulation)
- **Blue edges**: present only in the reference group
- **Orange edges**: present only in the case group
- **Red edges**: sign-switch edges (correlation reversed between groups)

**2. `triage/priority_rewired_nodes.csv`** — open in Excel, LibreOffice Calc, or any spreadsheet tool

This is the ranked list of metabolites most reorganised between the two states. The top of this list is your primary candidate set for hypothesis generation.

---

## How to read priority_rewired_nodes.csv

| Column | What it means |
|--------|---------------|
| `feature_id` | Metabolite identifier (same as in your feature table) |
| `rewiring_score` | Total number of edges that changed between groups: state-specific edges + sign-switch edges. Higher = more rewired. |
| `rewiring_fraction` | Fraction of this metabolite's connections that changed (0 = fully conserved neighbourhood, 1 = completely rewired neighbourhood) |
| `sign_switch_edges` | Number of edges that reversed correlation sign between groups (positive → negative or vice versa) |
| `shared_edges` | Number of edges conserved between both groups |
| `triage_score` | Combined score: 70% network centrality (topology_priority_score) + 30% phenotype association strength (1 − q_pheno). Range 0–1. Higher = more important and more strongly associated. |
| `q_pheno` | FDR-corrected p-value for the association between this feature's abundance and group identity |
| `log2FC` | log2(case mean / reference mean); positive = higher in case group |

Sort by `triage_score` (descending) to find features that are both strongly associated with the phenotype and topologically important in their network. Sort by `rewiring_score` (descending) to find features whose correlation context changed the most regardless of their centrality.

---

## How to read rewiring_summary.csv

This file contains one row with global statistics for the entire rewiring comparison.

| Column | What it means |
|--------|---------------|
| `n_shared_edges` | Edges conserved in both group networks |
| `n_state_A_only` | Edges present only in the reference group network |
| `n_state_B_only` | Edges present only in the case group network |
| `n_sign_switch_edges` | Edges where the correlation sign reversed between groups |
| `rewiring_proportion` | Fraction of all edges in the union network that are state-specific: (state_A_only + state_B_only) / total_union_edges |

A `rewiring_proportion` near 0 means the two networks are structurally very similar despite any abundance differences. A value above 0.5 means more than half of all co-regulation relationships changed between states — a signal of substantial network reorganization.

---

## How to read the state-specific network files

For each group (e.g., `networks/case/` and `networks/reference/`), PhenoRewire writes:

**`node_metrics.csv`** — one row per metabolite feature in that group's network

| Column | What it means |
|--------|---------------|
| `feature_id` | Metabolite identifier |
| `degree` | Number of co-regulation edges this metabolite has in this group |
| `degree_centrality` | Degree normalised by the maximum possible (range 0–1) |
| `weighted_degree_abs` | Sum of absolute Spearman rho across all edges of this node; measures total co-regulation strength |
| `betweenness_centrality` | How often this node sits on the shortest path between other nodes; high values indicate a hub or bridge |
| `eigenvector_centrality` | Importance weighted by the importance of neighbours; high = well-connected to other well-connected nodes |
| `community` | Louvain community membership ID (integer); nodes in the same community tend to co-vary together |
| `q_pheno` | FDR-corrected p-value for association with group identity (from selection step) |

**`edges.csv`** — one row per edge in that group's network

| Column | What it means |
|--------|---------------|
| `node_a`, `node_b` | The two metabolites connected by this edge |
| `rho` | Spearman correlation coefficient (positive = co-increase, negative = inverse relationship) |
| `q_edge` | FDR-corrected p-value for the correlation |

---

## Sign-switching edges — biological interpretation

A sign-switch edge is a co-regulation relationship that was positive in one state and negative in the other: for example, two metabolites that tend to increase together in the reference group but diverge (one increases while the other decreases) in the case group. This is qualitatively different from simply losing or gaining a correlation — it indicates a fundamental reorganization of the relationship between two metabolites.

Sign-switch edges are high-priority candidates for hypothesis generation because they suggest that the regulatory logic governing the relationship between two metabolites changed between biological states. Look for sign-switch edges involving metabolites with high `triage_score` or known pathway membership.

---

## What final_report.md contains

`final_report.md` is a plain-text narrative summary of the run. It includes:
- The config parameters used (FDR alpha, permutation count, network threshold)
- Feature selection statistics: how many features were tested, how many passed FDR, whether adaptive threshold relaxation was triggered
- Network statistics: number of nodes and edges per group, network modularity, Louvain stability
- Rewiring summary: counts of shared, state-specific, and sign-switch edges
- Top-ranked rewired features by triage_score

Use `final_report.md` for a quick orientation to any run and as the basis for the methods section of a paper. The companion `final_report.json` contains the same information in machine-readable format.

---

## Important — what PhenoRewire is NOT

- PhenoRewire identifies **correlational co-regulation patterns**, not causal relationships. A feature that appears in the top of `priority_rewired_nodes.csv` changed its correlation context between groups — this does not prove that it causes, mediates, or is responsible for the biological difference between groups.

- **Rewiring indicates reorganisation of context**, not effect size. A feature can be highly rewired but show no significant change in mean abundance, and vice versa. These two pieces of information are complementary, not redundant.

- **Use PhenoRewire outputs as hypotheses to test** with targeted experiments, mechanistic models, or orthogonal data (e.g., pathway enrichment, isotope labelling, enzyme activity assays). The pipeline provides a prioritised candidate list, not a causal explanation.
