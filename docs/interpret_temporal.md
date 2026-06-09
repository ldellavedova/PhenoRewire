# How to interpret PhenoRewire results — temporal mode

This guide is written for metabolomics researchers who may have no bioinformatics background. It explains what temporal mode does, which files to open first, and what each column means.

---

## What temporal mode does

Temporal mode compares two timepoints within a single phenotype — for example, 3h vs 6h within a control group, or baseline vs treatment-end within a patient cohort. Unlike phenotype mode, the question is not "which metabolites differ between two groups" but "which metabolites changed their co-regulation pattern as time progressed." Feature selection identifies metabolites whose abundance is associated with the time variable (treated as continuous). PhenoRewire then builds one correlation network per timepoint using only those selected features and compares the two networks edge by edge to quantify how the metabolic co-regulation structure evolved over time.

**Key difference from phenotype mode:** the association variable is the numeric time value, not a binary group label. Directionality is reported as increasing or decreasing abundance at the later timepoint, rather than as a log2 fold change between groups.

---

## Start here — your two primary output files

After every temporal run, open these two files first:

**1. `triage/rewiring_network_temporal.graphml`** — open in [Cytoscape](https://cytoscape.org) (free download)

Nodes are metabolite features. Node size is proportional to `rewiring_score` — larger nodes changed more between timepoints. Edge color distinguishes:
- **Grey edges**: shared between both timepoints (conserved co-regulation)
- **Blue edges**: present only at timepoint T1
- **Orange edges**: present only at timepoint T2
- **Red edges**: sign-switch edges (correlation reversed between timepoints)

**2. `triage/priority_rewired_nodes.csv`** — open in Excel, LibreOffice Calc, or any spreadsheet tool

This is the ranked list of metabolites most reorganised between the two timepoints. The top of this list is your primary candidate set for hypothesis generation about time-dependent metabolic rewiring.

---

## How to read priority_rewired_nodes.csv

| Column | What it means |
|--------|---------------|
| `feature_id` | Metabolite identifier (same as in your feature table) |
| `rewiring_score` | Total number of edges that changed between timepoints: timepoint-specific edges + sign-switch edges. Higher = more rewired. |
| `rewiring_fraction` | Fraction of this metabolite's connections that changed (0 = fully conserved neighbourhood, 1 = completely rewired neighbourhood) |
| `sign_switch_edges` | Number of edges that reversed correlation sign between timepoints |
| `shared_edges` | Number of edges conserved at both timepoints |
| `triage_score` | Combined score: 70% network centrality (topology_priority_score) + 30% time-association strength (1 − q_time). Range 0–1. |
| `q_time` | FDR-corrected p-value for the association between this feature's abundance and the time variable |
| `direction` | Whether this feature increased or decreased at the later timepoint (e.g., `higher_at_T2` or `lower_at_T2`) |

Sort by `triage_score` (descending) to find features that are both strongly associated with time and topologically important. Sort by `rewiring_score` (descending) to find features whose co-regulation context changed the most regardless of their centrality.

---

## How to read rewiring_summary.csv

This file contains one row with global statistics for the entire rewiring comparison.

| Column | What it means |
|--------|---------------|
| `n_shared_edges` | Edges conserved at both timepoints |
| `n_state_A_only` | Edges present only in the T1 network |
| `n_state_B_only` | Edges present only in the T2 network |
| `n_sign_switch_edges` | Edges where the correlation sign reversed between timepoints |
| `rewiring_proportion` | Fraction of all edges in the union network that are timepoint-specific: (T1_only + T2_only) / total_union_edges |

A `rewiring_proportion` near 0 means the co-regulation structure was stable over time. A high value indicates substantial temporal reorganization of the metabolic network.

---

## How to read the timepoint-specific network files

For each timepoint (e.g., `networks/T1/` and `networks/T2/`), PhenoRewire writes:

**`node_metrics.csv`** — one row per metabolite feature at that timepoint

| Column | What it means |
|--------|---------------|
| `feature_id` | Metabolite identifier |
| `degree` | Number of co-regulation edges at this timepoint |
| `degree_centrality` | Degree normalised by the maximum possible (range 0–1) |
| `weighted_degree_abs` | Sum of absolute Spearman rho across all edges; measures total co-regulation strength |
| `betweenness_centrality` | How often this node sits on the shortest path between others; high values indicate a hub |
| `eigenvector_centrality` | Importance weighted by the importance of neighbours |
| `community` | Louvain community membership ID; nodes in the same community tend to co-vary together |
| `q_time` | FDR-corrected p-value for association with the time variable |
| `direction` | `higher_at_T2` or `lower_at_T2` |

**`edges.csv`** — one row per edge in that timepoint's network

| Column | What it means |
|--------|---------------|
| `node_a`, `node_b` | The two metabolites connected by this edge |
| `rho` | Spearman correlation coefficient (positive = co-increase, negative = inverse relationship) |
| `q_edge` | FDR-corrected p-value for the correlation |

---

## Sign-switching edges — biological interpretation

A sign-switch edge indicates that two metabolites co-increased at one timepoint but diverged (one increased, the other decreased) at the other timepoint, or vice versa. In a temporal context, this is particularly informative: it suggests that the metabolic relationship between two features reversed as time progressed, which may reflect a shift in the regulatory mechanism governing both metabolites.

Look for sign-switch edges involving metabolites with high `triage_score` or known involvement in time-sensitive processes (e.g., circadian metabolism, drug biotransformation, substrate-product relationships).

---

## What final_report.md contains

`final_report.md` is a plain-text narrative summary of the run. It includes:
- The config parameters used
- Feature selection statistics: how many features were tested, how many passed FDR, whether adaptive threshold relaxation was triggered
- Network statistics per timepoint: number of nodes, edges, modularity, Louvain stability
- Rewiring summary: counts of shared, timepoint-specific, and sign-switch edges
- Top-ranked rewired features by triage_score

Use `final_report.md` for a quick orientation to any run and as the basis for the methods section of a paper. The companion `final_report.json` contains the same information in machine-readable format.

---

## Important — what PhenoRewire is NOT

- PhenoRewire identifies **correlational co-regulation patterns over time**, not causal dynamics. A metabolite at the top of `priority_rewired_nodes.csv` changed its correlation context between timepoints — this does not prove that it drives, mediates, or is responsible for the temporal progression.

- **Temporal rewiring is not the same as a time-course abundance change.** A feature can show strong rewiring with minimal change in mean abundance, and a feature can change dramatically in abundance while maintaining the same co-regulation partners. These are independent pieces of information.

- **Use PhenoRewire outputs as hypotheses to test** with targeted experiments, isotope tracing, or mechanistic models. The pipeline provides a prioritised candidate list for time-dependent metabolic reorganisation, not a causal timeline.
