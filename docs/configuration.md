# Configuration reference

PhenoRewire is driven by one YAML file. It is validated before any data is
read, and an invalid key or value stops the run with a message naming it.
Relative paths are resolved against the directory that holds the YAML file.
`phenorewire --config my_config.yaml --check` validates the config and reads
the data it points to without running the analysis.

This page lists every key the validator accepts. The two annotated templates
at the repository root, `config_phenotype.yaml` and `config_temporal.yaml`,
are complete working examples.

---

## A minimal phenotype config

Only the keys without a default are required. This config compares two
groups and leaves everything else at its default:

```yaml
DATA_MATRIX: "data/feature_table.csv"
METADATA: "data/metadata.csv"
OUTDIR: "results/minimal"
ANALYSIS_MODE: "phenotype"

GROUP_DEFINITION:
  Reference:
    group_label: "reference"
  Case:
    group_label: "case"
PHENO_GROUP_REF: "Reference"
PHENO_GROUP_CASE: "Case"

FEATURE_ID_COL: "feature_id"
MZ_COL: "mz"
RT_COL: "rt"
INTENSITY_REGEX: "^Sample_"
META_SAMPLE_COL: "sample_id"
```

The test suite validates this block against the bundled example data, so it
cannot drift out of date.

---

## Three keys that need explaining

### `INTENSITY_REGEX`

A regular expression that picks the intensity columns out of the feature
table. Every column it matches is treated as one sample; every other column
is treated as annotation.

The sample ID is derived from the column name: when the match sits at the
end of the name it is removed as a suffix (as are the MZmine suffixes `Peak
area`, `Peak height`, `Area` and `Height`); a match elsewhere leaves the name
as it is. File extensions such as `.mzML` are dropped from both the column
names and the metadata IDs before they are compared.

| Feature-table columns | `INTENSITY_REGEX` | Sample IDs | Metadata values that match |
|---|---|---|---|
| `Sample_01`, `Sample_02` | `"^Sample_"` | `Sample_01`, `Sample_02` | `Sample_01` |
| `S1.mzML Peak area` (MZmine) | `"Peak area$"` | `S1` | `S1` or `S1.mzML` |
| `S1 Area` | `"Area$"` | `S1` | `S1` |

If no column matches, the run stops with "No intensity columns found".

### `GROUP_DEFINITION`

A mapping from a group name you choose to the metadata values that define
it. A sample belongs to a group when every condition holds:

```yaml
GROUP_DEFINITION:
  Full_community_3h:          # group name used in outputs
    all_strains: 1            # metadata column: required value
    timepoint: "3h"
  Without_W_anitrata_3h:
    no_anitrata: 1
    timepoint: "3h"
```

Values are compared tolerantly: `1`, `"1"` and `1.0` match each other, as do `true` and
`"yes"`. A sample matching several groups goes to the first one listed;
samples that match no group are ignored. `PHENO_GROUP_REF` and
`PHENO_GROUP_CASE` name the two groups to contrast; the reference is coded 0
and the case 1, so a positive correlation or log2 fold change means higher in
the case. When more than two groups match samples, these two keys select the
contrast.

### `ANALYSIS_MODE`

| Value | What runs |
|---|---|
| `phenotype` | Two groups compared |
| `temporal` | Two timepoints compared within one group (needs `TEMPORAL_CORRELATION.enabled: true`) |
| `both` | Both of the above in one run |
| `auto` | Whatever the config and data support |

A mode that was asked for and cannot run stops the run; it is never skipped
in silence.

---

## All keys

### Inputs and outputs

| Key | Default | Meaning |
|---|---|---|
| `DATA_MATRIX` | required | Feature table, one row per feature |
| `METADATA` | required | One row per sample |
| `OUTDIR` | required | Output directory; created when the run starts |
| `ANALYSIS_MODE` | `auto` | See above |

### Groups and timepoints

| Key | Default | Meaning |
|---|---|---|
| `GROUP_DEFINITION` | required | See above |
| `PHENO_GROUP_REF`, `PHENO_GROUP_CASE` | none | The two groups to contrast; set both or neither |
| `TIMEPOINT_T1`, `TIMEPOINT_T2` | 3, 6 | Timepoints compared in temporal mode; T1 < T2 |

### Column names

| Key | Default | Meaning |
|---|---|---|
| `FEATURE_ID_COL` | required | Unique feature identifier |
| `MZ_COL`, `RT_COL` | required | m/z and retention time |
| `NAME_COL` | none | Optional compound name or annotation |
| `INTENSITY_REGEX` | required | See above |
| `META_SAMPLE_COL` | required | Sample identifier in the metadata |
| `META_TIME_COL` | none | Time column; fallback for `TEMPORAL_CORRELATION.time_column` |

### Preprocessing

| Key | Default | Meaning |
|---|---|---|
| `MIN_NONZERO_INTENSITY` | 0.0 | Values at or below this count as absent |
| `MIN_PRESENCE_GLOBAL` | 0.25 | Keep features present in at least this fraction of samples |
| `NORM_METHOD` | `median` | Sample normalisation: `median`, `pqn` (probabilistic quotient), `tic` (total signal) or `none` (data already normalised). Which to choose: [`prepare_your_data.md`](prepare_your_data.md#choosing-the-normalisation-norm_method) |
| `LOG_TRANSFORM` | true | Apply log2(x + 1) after normalisation |

Features with the same intensity in every sample (typically values imputed
at a detection limit) are removed: they carry no sample-to-sample
information and would otherwise correlate perfectly with each other.

### Connectivity-based selection

Runs automatically, alongside abundance-based selection, whenever both states
of a comparison have at least `CONNECTIVITY.min_n_per_state` samples. Below
that, features are selected by abundance only and the log and report say why.
The settings rarely need changing.

| Key | Default | Meaning |
|---|---|---|
| `CONNECTIVITY.top_k` | 10 | A feature's score is the mean of its `top_k` largest changes in Spearman rho to any partner |
| `CONNECTIVITY.n_permutations` | 200 | Label permutations for its null |
| `CONNECTIVITY.fdr_alpha` | 0.05 | BH threshold |
| `CONNECTIVITY.min_n_per_state` | 14 | Samples needed in each state; minimum allowed value 10 |

When it runs, every selected feature carries `selection_route`
(`abundance`, `connectivity` or `both`). Features with route `connectivity`
changed their co-variation partners without a detectable change in abundance:
a differential-abundance analysis cannot report them.
`connectivity_features_all.csv` in the selection folder holds every feature's
score, p-values and q-value.

### Feature selection (stage 1)

| Key | Default | Meaning |
|---|---|---|
| `FDR_ALPHA` | 0.05 | BH threshold for selecting features |
| `N_PERMUTATIONS_PHENO` | 1000 | Permutations per feature. The smallest attainable p is 1/(n+1); below ~1000 it approaches the BH threshold and suppresses selection |
| `ADAPTIVE_SELECTION_THRESHOLDS` | false | If too few features pass `FDR_ALPHA`, relax it stepwise up to the ceiling. Flagged prominently in the report when used |
| `SELECTION_MIN_FEATURES` | 5 | Target the relaxation aims for |
| `SELECTION_FDR_ALPHA_CEILING` | 0.10 | Highest alpha the relaxation may reach |
| `MIN_FEATURES_HARD_STOP` | 5 | Stop if fewer features than this are selected |
| `RANDOM_SEED` | 42 | Seed for every permutation in the run. Selection permutations are drawn from one generator consumed feature by feature, the scheme all published PhenoRewire analyses used, so a published analysis is reproduced exactly with its seed |

### Correlation networks (stage 2)

| Key | Default | Meaning |
|---|---|---|
| `CORR_ABS_THRESHOLD` | 0.6 | Minimum \|Spearman rho\| for an edge |
| `CORR_FDR_ALPHA` | 0.05 | BH threshold across all pairs in one state |
| `ADAPTIVE_NETWORK_THRESHOLDS` | true | Search a small grid around the two thresholds above for a network that is neither empty nor saturated; the thresholds applied are recorded in `network_summary.csv` |
| `NETWORK_MIN_EDGES` | 3 | Lower edge target for the adaptive search |
| `NETWORK_MAX_DENSITY` | 0.25 | Upper density target for the adaptive search |
| `NETWORK_ABS_THRESHOLD_FLOOR` | 0.5 | Lowest \|rho\| the search may use |
| `NETWORK_FDR_ALPHA_CEILING` | 0.2 | Highest FDR the search may use |
| `NETWORK_MIN_EDGES_HARD_STOP` | 10 | Stop if a final network has fewer edges |
| `NETWORK_MIN_NODES_HARD_STOP` | 5 | Stop if a final network has fewer nodes |
| `LOUVAIN_STABILITY_CHECK` | false | Run Louvain ten times and report agreement (needs scikit-learn) |

When the adaptive search applies different thresholds to the two states of
a comparison, the rewiring null uses each state's own thresholds, so the
comparison against chance is unaffected by the asymmetry.

### Rewiring (stage 3) and null models

| Key | Default | Meaning |
|---|---|---|
| `SIGN_SWITCH_MIN_R` | 0.5 | A shared edge counts as a sign switch only if \|rho\| >= this in both states |
| `NULL_MODEL.enabled` | true | Compute the null models below. They add files and never change any other output |
| `NULL_MODEL.n_permutations` | 200 | Replicates for the edge null of each network and the label-permutation null of each comparison (minimum 20) |
| `NULL_MODEL.n_configuration_models` | 50 | Degree-preserving rewirings used to put modularity in context |

### Ranking

| Key | Default | Meaning |
|---|---|---|
| `REPORT_TOP_N` | 15 | Metabolites listed in `REPORT.md` and `report.json`; the CSV tables always list all of them |
| `TRIAGE_TOPOLOGY_WEIGHT` | 0.7 | Weight of network position in `triage_score` |
| `TRIAGE_SELECTION_WEIGHT` | 0.3 | Weight of selection strength; the two must sum to 1 |

### Temporal mode

| Key | Default | Meaning |
|---|---|---|
| `TEMPORAL_CORRELATION.enabled` | false | Turn temporal mode on |
| `TEMPORAL_CORRELATION.phenotype` | none | The group whose time course is analysed |
| `TEMPORAL_CORRELATION.time_column` | none | Metadata column with times; "3h" and "3" both parse as 3 |
| `TEMPORAL_CORRELATION.n_permutations` | 1000 | As `N_PERMUTATIONS_PHENO` |
| `TEMPORAL_CORRELATION.fdr_alpha` | 0.05 | As `FDR_ALPHA` |
| `TEMPORAL_CORRELATION.min_samples` | 3 | Skip with a warning below this many samples in the group |

Temporal selection correlates each feature with time across all timepoints
of the group; the networks are then built at `TIMEPOINT_T1` and
`TIMEPOINT_T2`.

### Molecular-network annotation (`--run-2dnetwork`)

| Key | Default | Meaning |
|---|---|---|
| `SPECREBOOT_NETWORK` | none | GraphML molecular network whose node IDs match `FEATURE_ID_COL` |
| `PHENOREWIRE_NETWORK` | none | Rewiring GraphML to start from; defaults to this run's |
| `MN_HOP_DEPTH` | 2 | Neighbourhood depth searched |
| `MN_MIN_COSINE` | 0.70 | Minimum spectral similarity of traversed edges |
| `ANNOTATION_TOP_N` | 20 | Prioritised features to expand |

---

## Sample size

There is no setting for it, but it matters more than any setting here. Each
run reports where the smaller state of every comparison sits against the
measured operating range: the input check prints it at the start of the run,
and `REPORT.md` and `report.json` (`operating_range`) repeat it. See
[`prepare_your_data.md`](prepare_your_data.md#sample-size) for what each range
means, and the README for the benchmark behind the numbers.
