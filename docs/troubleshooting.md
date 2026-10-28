# Troubleshooting

Start with the input check: it reads your tables exactly as a run would and names what is wrong.

```bash
phenorewire --config my_config.yaml --check
```

When a run stops, the message says why and how to fix it. The full log, with timestamps, is in `OUTDIR/run_info/run.log`; `--log-level DEBUG` also prints the Python traceback. Attach `run.log` when you report a problem.

---

## Installation

### `phenorewire` command not found

Activate the environment you installed into, then install again:

```bash
conda activate phenorewire   # or: source .venv/bin/activate
pip install -e .
```

### Python version

PhenoRewire needs Python ≥ 3.9 and is tested on 3.9–3.12. `environment.yml` pins 3.11.

### `--plot-networks` or the figure script fails with ImportError

```bash
pip install "phenorewire[viz]"
```

---

## The config

### `ERROR: the config file ... has a problem`

Each line names the key and what is wrong with it:

| Message | Fix |
|---|---|
| `DATA_MATRIX: file not found: ...` | The path is wrong. Relative paths are read from the folder that holds the config file, not from where you run the command |
| `LOG_TRANSFORMS: Extra inputs are not permitted` | A misspelled key. Every valid key is in [`configuration.md`](configuration.md) |
| `TIMEPOINT_T1 (6.0) must be less than TIMEPOINT_T2 (3.0)` | Swap the two values |

---

## The input data

These stop the run at step 1, before any statistics.

### No sample of the feature table has a row in the metadata

The message shows sample names from both tables. Common causes:

- `INTENSITY_REGEX` leaves part of the column name in the sample name (`Sample_01 Peak area` must become `Sample_01`); see [`configuration.md`](configuration.md#intensity_regex).
- `META_SAMPLE_COL` points at the wrong column.
- Different spelling, for example `S01` vs `Sample_01`.

### Some samples of the feature table have no row in the metadata

Every intensity column needs a metadata row. Add the rows, or remove those columns (blanks, QCs) from the feature table, or narrow `INTENSITY_REGEX` so it does not select them.

### The intensity values are negative

The table has been autoscaled, mean-centred or log-transformed. Provide the intensities before scaling: see [`prepare_your_data.md`](prepare_your_data.md#intensities-raw-non-negative-not-scaled).

### Features hold text instead of intensities

Usually binary variables (`TRUE`/`FALSE`, `VRAI`/`FAUX`). They are removed; remove them from the table to silence the warning. See [`prepare_your_data.md`](prepare_your_data.md#binary-variables).

### The table looks already log-transformed

Set `LOG_TRANSFORM: false`, or provide the values before the log transform.

### No sample matches any group / GROUP_DEFINITION refers to columns that do not exist

Values in `GROUP_DEFINITION` must match the metadata exactly, including case (`case` is not `Case`) and type (the YAML `1` matches a metadata `1`). The check lists the metadata columns it found.

### N samples match no group in GROUP_DEFINITION and are not used

A warning, not an error, and often correct: it is what keeps QC injections, blanks and the conditions a given comparison leaves out from entering the networks. Read the names it lists. If a sample you meant to analyse is among them, a value in `GROUP_DEFINITION` does not match the metadata — check case and type, as above. If the count is larger than you expect, the groups are matching fewer rows than you think, and the comparison is running on fewer samples than the design has.

Blanks and QC samples are better removed from the feature table altogether, because the presence filter and the normalisation reference are computed before groups are assigned: see [`prepare_your_data.md`](prepare_your_data.md#quality-control-samples-blanks-and-pools).

### TIMEPOINT_T2 = 10 has no samples

The check lists the timepoints that exist in the group. Set `TIMEPOINT_T1` and `TIMEPOINT_T2` to two of them.

---

## During the run

### Too few selected features

**Message:** `only N feature(s) selected, fewer than MIN_FEATURES_HARD_STOP`.

Few metabolites differ between the states at the chosen FDR. Options: more samples per state; `FDR_ALPHA: 0.1`; `ADAPTIVE_SELECTION_THRESHOLDS: true` (the relaxation is reported); for a small exploratory dataset, `MIN_FEATURES_HARD_STOP: 2`.

### Network has only N edges

**Message:** `Network 'X' has only N edge(s) after all threshold adjustments`.

The default `CORR_ABS_THRESHOLD` (0.6) can be too strict for large, heterogeneous cohorts, where real correlations between metabolites are weaker. Lower `CORR_ABS_THRESHOLD` and `NETWORK_ABS_THRESHOLD_FLOOR` (for example to 0.3), or raise `CORR_FDR_ALPHA`. For small demo datasets, `NETWORK_MIN_EDGES_HARD_STOP: 1`.

### The rewiring is not distinguishable from chance

Not an error: the two networks differ no more than two arbitrary halves of the same samples would. Below 14 samples per state this is the expected result even when something has changed; see the sample-size table in [`reading_the_results.md`](reading_the_results.md#sample-size).

---

## Molecular-network annotation (`--run-2dnetwork`)

### Zero neighbours, or many `in_specreboot = False` rows

The molecular-network GraphML must use the same node IDs as the `feature_id` of your feature table.

### Propagation goes through every edge regardless of similarity

The molecular-network GraphML needs a numeric similarity attribute on its edges, such as `cosine`.

---

## Verify the installation

```bash
phenorewire --help
python -m pytest -q dev/tests/
```
