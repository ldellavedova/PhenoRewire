# PhenoRewire — Example Walkthrough

This directory contains a self-contained example to get you started and to explain how to interpret PhenoRewire outputs.

## Dataset description

| | |
| --- | --- |
| Features | 22 metabolites (20 biologically patterned + 2 noise) |
| Samples | 20 (10 reference, 10 case) |
| Phenotype | two-group comparison (reference vs case) |

The synthetic data embeds a specific biological scenario:

- **Nucleotide energy module** (ATP, ADP, AMP, adenosine): strongly co-regulated in the reference group; this correlation structure dissolves in the case group
- **Case-specific feature**: one feature gains a new correlation partner only in the case group
- Several features are differentially abundant between groups

## Step 1 — Run the pipeline

From the repository root:

```bash
phenorewire --config dev/benchmarks/fixtures/config_example.yaml
```

Results are written to `results/example/`. The run takes about 30–60 seconds.

## Step 2 — Read the summary first

Open `results/example/report/final_report.md`.

The **Executive Summary** section gives you the headline numbers immediately:

```
Phenotype comparison retained 9 phenotype-associated features.
Phenotype rewiring summary: shared=0, condition_specific=7, sign_switch=0.
```

Interpretation:
- 9 of 22 features were significantly correlated with the phenotype after FDR correction (q ≤ 0.10)
- Of the 7 total edges across both networks, **0 are shared** between reference and case — the entire correlation structure is rewired
- No sign switches in this example (an edge that was positive in one group is not negative in the other)

The **rewiring_proportion = 1.0** in `rewiring_summary.csv` confirms this: 100% of all connections changed between states.

## Step 3 — Identify the most rewired metabolites

Open `results/example/triage/priority_rewired_nodes.csv`.

The top rows will look like this:

| feature_id | name | rewiring_score | rewiring_fraction | Reference_only_edges | Case_only_edges | sign_switch_edges |
| --- | --- | --- | --- | --- | --- | --- |
| F007 | ADP | 4 | 1.0 | 3 | 1 | 0 |
| F005 | adenosine | 3 | 1.0 | 3 | 0 | 0 |
| F006 | AMP | 3 | 1.0 | 3 | 0 | 0 |
| F008 | ATP | 3 | 1.0 | 3 | 0 | 0 |

Reading this table:

- **ADP (F007)** has the highest `rewiring_score` (4): it had 3 edges only in the reference network and 1 edge only in the case network. None of its connections are shared between groups.
- **adenosine, AMP, ATP** each lost 3 reference-specific connections and gained none in the case group — the nucleotide energy module essentially disappeared in the case condition.
- `rewiring_fraction = 1.0` for all top nodes means every single one of their connections changed state. These are not merely losing one edge; they are completely restructured.

> **Biological hypothesis**: the loss of co-regulation among ATP, ADP, AMP and adenosine in the case group suggests disruption of energy charge homeostasis. These four metabolites normally track together through adenylate kinase and ATPases; when they decouple, it signals impaired energy buffering.

## Step 4 — Compare the two state-specific networks

Open the per-state networks in Cytoscape (or any GraphML viewer):

- `results/example/network_Reference/network.graphml` — 9 nodes, 6 edges
- `results/example/network_Case/network.graphml` — 9 nodes, 1 edge

**Reference network**: ATP, ADP, AMP, and adenosine form a tight, fully connected module (community 1). The 6 edges reflect their co-regulation under normal energy metabolism.

**Case network**: the nucleotide module has collapsed. Only one edge remains (ADP ↔ F021). The network is nearly star-free.

This collapse is the core finding: the correlation structure that defines energy metabolism in the reference group is absent in the case group.

## Step 5 — Check individual feature statistics

Open `results/example/phenotype_selection/selected_pheno_features_with_direction.csv`.

Key columns:

| Column | Meaning |
| --- | --- |
| `r_pheno` | Spearman correlation with phenotype assignment (positive = higher in case) |
| `q_pheno` | BH-corrected FDR q-value from permutation test |
| `log2FC_case_vs_ref` | Mean(case) − Mean(ref) in log-scale; positive = up in case |
| `direction` | `up_in_Case` or `up_in_Reference` |

Example rows:

| name | r_pheno | q_pheno | log2FC_case_vs_ref | direction |
| --- | --- | --- | --- | --- |
| phenylalanine | +0.87 | 0.012 | +2.70 | up_in_Case |
| adenosine | −0.76 | 0.012 | −1.53 | up_in_Reference |
| ATP | −0.76 | 0.012 | −1.22 | up_in_Reference |

Phenylalanine (F013) is strongly elevated in the case group. The nucleotide metabolites are elevated in the reference group.

## Step 6 — Visualise in Cytoscape (optional)

1. Open Cytoscape → File → Import → Network from File → select `results/example/triage/rewiring_network_phenotype.graphml`
2. Style → Edge colour mapped to `state`:
   - `Reference_only` → blue
   - `Case_only` → red
   - `shared` → grey
3. Style → Node size mapped to `rewiring_score`

The result will immediately show which metabolites are hubs of change.

## Summary of output files

```
results/example/
├── triage/                          ← START HERE
│   ├── final_report.md              ← narrative summary (this walkthrough refers to it)
│   ├── priority_rewired_nodes.csv   ← top rewired metabolites
│   ├── rewiring_summary.csv         ← one-row global stats
│   └── rewiring_network_phenotype.graphml  ← open in Cytoscape
│
├── phenotype_selection/
│   └── selected_pheno_features_with_direction.csv   ← per-feature stats + direction
│
├── network_Reference/               ← reference-group network (detail)
├── network_Case/                    ← case-group network (detail)
└── rewiring_pheno/                  ← full rewiring tables (detail)
```

Files marked **detail** are for advanced users who want to inspect individual edges, community membership, or full node-metric tables.
