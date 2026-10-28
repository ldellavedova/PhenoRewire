# Reading the results

This guide explains what PhenoRewire writes and what every column means. It is the same for a comparison of two groups (phenotype mode) and of two timepoints within one group (temporal mode); where the two differ, both are given.

No bioinformatics background is assumed. For a worked example on the bundled data, see [`tutorials/walkthrough.md`](tutorials/walkthrough.md).

---

## The output folder

```text
OUTDIR/
├── REPORT.md                        start here
├── report.json                      everything in REPORT.md, machine-readable
├── <comparison>/                    one folder per comparison
│   ├── rewiring_vs_chance.csv       1. do the two networks differ more than chance?
│   ├── rewired_metabolites.csv      2. one row per metabolite, strongest evidence first
│   ├── selected_metabolites.csv     3. which metabolites entered the networks, and why
│   ├── rewiring_network.graphml     4. both networks in one file, for Cytoscape
│   └── details/
│       ├── selection/               the statistics of every feature tested
│       ├── network_<state>/         one folder per network
│       ├── rewiring/                every correlation with its state
│       └── samples_used.csv         the samples in each state
└── run_info/                        log, input check, configuration, preprocessed data
```

The comparison folder is named after the two states: `Reference_vs_Case/` in phenotype mode, `<group>_<T1>_vs_<T2>/` in temporal mode (for example `Treated_0h_vs_10h/`). A run with `ANALYSIS_MODE: both` writes one of each.

**Read the four files at the top of the comparison folder, in order.** Everything in `details/` and `run_info/` is there to check, reuse or plot the results; you do not need it to interpret them.

---

## REPORT.md

A one-page summary of the run, meant to be read first:

- **Result** — one sentence per comparison: do the two networks differ more than chance, and how many metabolites rewire more than chance.
- **Warnings** — anything that limits the interpretation: a sample size below the operating range, a rewiring not distinguishable from chance, a relaxed FDR threshold, problems found in the input data.
- **Start here** — the files to open, and a table of the most rewired metabolites.
- **Details** — how many metabolites were selected and by which route, the sample size per state, and one line per network (metabolites, correlations, thresholds, comparison with chance).
- **How the numbers are computed** and **Reproducibility** — the formulas, the seed and every threshold applied. Use them for the methods section of a paper.

---

## 1. `rewiring_vs_chance.csv` — is the rewiring real?

Read this before any ranking. Two networks built from different samples always differ somewhat; the question is whether they differ more than two arbitrary halves of the same samples would.

PhenoRewire answers it by permutation: the state labels are shuffled across the pooled samples (each metabolite first centred within its own state, so that abundance differences cannot create correlations), both networks are rebuilt with the same thresholds, and the state-specific correlations are counted. One row:

| Column | Meaning |
|---|---|
| `observed_state_specific_edges` | Correlations present in exactly one of the two networks |
| `null_mean`, `null_sd`, `null_p95` | The same count over the shuffled splits: its mean, standard deviation and 95th percentile |
| `obs_over_null` | Observed divided by `null_mean`; about 1 means no more than chance |
| `p_empirical` | (number of shuffles with a count ≥ observed + 1) / (shuffles + 1) |
| `n_a`, `n_b`, `n_features` | Samples per state and metabolites in the networks |
| `applied_abs_r_thr_*`, `applied_fdr_alpha_*` | The thresholds of each network, reused for every shuffle |
| `n_permutations`, `centred_within_state` | How the null was built |

**If `p_empirical` > 0.05, stop here:** the rest describes these two networks, but not a difference between the states. With 200 permutations (the default) the smallest possible p is 0.005.

This file is written when `NULL_MODEL.enabled` is true (the default).

---

## 2. `rewired_metabolites.csv` — which metabolites rewired

One row per metabolite that entered the networks, most rewired first: sorted by `q_rewiring`, then by `rewiring_score`.

| Column | Meaning |
|---|---|
| `rank` | Position in this table |
| `feature_id`, `name` | The metabolite, as in your feature table |
| `rewiring_score` | Correlations present in one network only, plus correlations that changed sign. Higher = more rewired |
| `rewiring_fraction` | Share of this metabolite's correlations present in one network only (0 = all kept, 1 = all changed) |
| `q_rewiring` | Is this metabolite's number of state-specific correlations more than chance? Benjamini–Hochberg q across metabolites, from the same permutations as `rewiring_vs_chance.csv` |
| `more_than_chance` | `q_rewiring` ≤ 0.05 |
| `selection_route` | How it was selected: `abundance`, `connectivity` (its co-variation partners changed, its abundance did not: a differential-abundance analysis cannot report it) or `both`. Present when every state has at least 14 samples |
| `direction` | Which state has the higher mean abundance: `up_in_<group>` (phenotype) or `up_at_<timepoint>` (temporal) |
| `<A>_only_edges`, `<B>_only_edges` | Correlations present only in the network of state A, or of state B |
| `shared_edges` | Correlations present in both networks |
| `sign_switch_edges` | Correlations present in both networks with opposite sign (positive in one, negative in the other) |
| `expected_by_chance`, `p_rewiring` | The metabolite's mean state-specific count over the shuffles, and its empirical p |

Metabolites with `more_than_chance = True` are the list to take forward. The raw `rewiring_score` grows with the number of correlations a metabolite has; use `rewiring_fraction` to compare metabolites with very different numbers of partners.

---

## 3. `selected_metabolites.csv` — what entered the networks

Only the metabolites listed here are in the networks. A metabolite is selected when its abundance differs between the states, or — with at least 14 samples per state — when its co-variation partners differ even though its abundance does not.

| Column | Meaning |
|---|---|
| `feature_id`, `name`, `mz`, `rt` | The metabolite |
| `selection_route` | `abundance`, `connectivity` or `both` (see above) |
| `direction` / `time_direction` | Which state has the higher mean abundance |
| `mean_ref`, `mean_case`, `log2FC_case_vs_ref` | Phenotype mode: mean preprocessed (log2) abundance per group, and their difference, a log2 fold change. Positive = higher in the case group |
| `mean_t1`, `mean_t2`, `delta_t2_minus_t1` | Temporal mode: the same at the two timepoints compared |
| `r_pheno`, `q_pheno` / `r_time`, `q_time` | Spearman correlation of the abundance with the group label (phenotype) or with time across all timepoints (temporal), and its Benjamini–Hochberg q |
| `connectivity_score`, `q_connectivity` | The mean of the metabolite's 10 largest changes in correlation to any partner, and its q against shuffled labels |
| `q_selection` | The q-value of the route that selected the metabolite (the smaller one when both did) |

In temporal mode, selection uses every timepoint of the group, while the networks compare only `TIMEPOINT_T1` and `TIMEPOINT_T2`.

---

## 4. `rewiring_network.graphml` — the picture

Both networks in one file, for [Cytoscape](https://cytoscape.org) (File → Import → Network from File). Nodes are metabolites, with every column of `rewired_metabolites.csv` and `selected_metabolites.csv` as attributes. Edges are correlations, with:

| Edge attribute | Meaning |
|---|---|
| `state` | `shared`, `<A>_only` or `<B>_only` |
| `sign_switch` | True when the correlation changed sign |
| `r_<A>`, `r_<B>`, `q_<A>`, `q_<B>` | The correlation and its q in each network (empty where absent) |

A useful style: colour edges by `state`, make sign-switched edges dashed, size nodes by `rewiring_score`, and highlight nodes with `more_than_chance`.

**Sign switches.** A correlation that is positive in one state and negative in the other (two metabolites that rise together in one state and diverge in the other) is a stronger change than a correlation that appears or disappears. These pairs are often worth following up first.

---

## details/

### `selection/`

| File | Contents |
|---|---|
| `all_features_tested.csv` | Every feature after preprocessing, with its association with the state or with time: `r_pheno`/`r_time`, `p_empirical`, `q_pheno`/`q_time`, `significant` (selected by abundance), and the group means |
| `all_features_connectivity.csv` | Every feature's connectivity statistics: `connectivity_score`, its null mean and z-score, `p_connectivity`, `q_connectivity`, `selected_connectivity`. Written when every state has at least 14 samples |
| `selection_summary.json` | Counts, thresholds and whether the FDR threshold was relaxed adaptively |

### `network_<state>/`

| File | Contents |
|---|---|
| `nodes.csv` | One row per metabolite in this network, best first by `triage_score`, with `degree`, `degree_centrality`, `weighted_degree_abs` (sum of \|r\| over its correlations), `betweenness_centrality` (how often it lies between other metabolites), `eigenvector_centrality` (connected to well-connected metabolites), `community_id` (Louvain module), and the selection columns |
| `edges.csv` | One row per correlation: `node_a`, `node_b`, `r` (Spearman), `q` (Benjamini–Hochberg), `sign` |
| `communities.csv` | One row per Louvain module: size, internal correlations, members |
| `network_summary.csv` | Nodes, edges, density, modularity, and the thresholds applied |
| `null_model.csv` | This network against chance: correlations expected if the metabolites varied independently (`null_mean`, `p_empirical`), and `delta_modularity`, the modularity beyond what the degree sequence alone gives |
| `network.graphml` | This network alone, for Cytoscape |

`triage_score` = 0.7 × topology rank + 0.3 × (1 − q of selection): high for metabolites that are both central in this network and strongly selected. The weights are `TRIAGE_TOPOLOGY_WEIGHT` and `TRIAGE_SELECTION_WEIGHT`.

Raw `modularity` rises as networks get sparser; compare networks of different density by `delta_modularity`, not by `modularity`.

### `rewiring/`

| File | Contents |
|---|---|
| `edges.csv` | Every correlation present in at least one network: its r, q and sign in each state, its `state` and `sign_switch` |
| `rewiring_summary.csv` | One row: `n_shared_edges`, `n_<A>_only_edges`, `n_<B>_only_edges`, `n_sign_switch_edges`, `rewiring_proportion` (share of all correlations that are state-specific) and `sign_switch_proportion` |

### `samples_used.csv`

The metadata of the samples in the two networks, with their `state`.

---

## run_info/

| File | Contents |
|---|---|
| `run.log` | Everything the run did, with timestamps. Attach it when you report a problem |
| `data_check.txt` | The input check: samples matched, groups, timepoints, sample sizes, warnings |
| `config_used.yaml` | The configuration as applied, every default filled in. Rerun with it to reproduce the results |
| `samples.csv` | Every sample of the feature table with its group |
| `feature_annotation.csv` | ID, m/z, RT and name of the features kept after preprocessing |
| `preprocessed_matrix.csv` | The data the statistics ran on: normalised, log-transformed intensities, features × samples |

---

## Sample size

Every comparison is classified by the number of samples in its smaller state. The input check prints it at the start of the run, and `REPORT.md` repeats it.

| Samples in the smaller state | Status | What to expect |
|---|---|---|
| ≥ 14 | within the operating range | Metabolites are also selected by change in co-variation partners; real rewiring is usually detected |
| 10–13 | reduced | Selection by abundance only. What is called is generally reliable; much real rewiring is missed |
| < 10 | below | Networks are dominated by near-perfect rank correlations. Treat any result as preliminary |

An absence of rewiring below 14 samples per state is not evidence that nothing rewired.

---

## What PhenoRewire does not tell you

- **Correlation, not cause.** A rewired metabolite changed its co-variation partners between the states. That does not show that it causes, mediates or responds to the difference between them.
- **Rewiring is not effect size.** Two metabolites with the same fold change can be fully rewired or fully conserved, and a metabolite selected by connectivity can rewire without any change in abundance.
- **Over time, not dynamics.** In temporal mode the two networks are snapshots at two timepoints; PhenoRewire does not model the trajectory between them.
- **Hypotheses to test.** Follow the top metabolites up with targeted experiments, pathway knowledge, isotope labelling or orthogonal data.
