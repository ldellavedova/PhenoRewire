# Changelog

All notable changes to PhenoRewire are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [0.2.0] — unreleased

Every network and every comparison is now reported against chance; metabolites
whose co-variation partners change without a change in abundance are selected
automatically when the data allow; and the tool states the sample sizes at which
its results can be trusted. Earlier analyses are kept: abundance-based selection
reproduces every published analysis exactly, connectivity-based selection only
adds metabolites, and the null models only add files and report sections.

The tool is also easier to start with: every run checks the input data first
and stops with a readable message on data it cannot analyse, the terminal shows
numbered steps and ends by saying which file to open, and the output is
reorganised around four files per comparison and a one-page `REPORT.md`
(see **Output layout** under Changed).

### Why

A ground-truth benchmark showed that the raw count of state-specific edges is
not a reliable signal on its own. When a block of metabolites shifts in
abundance but keeps its correlation structure, stage 1 selects the block, and
some of its correlations then cross the edge threshold in one state and not
the other by chance alone. The raw count therefore reports rewiring that is
not there, and more of it as sample size grows. The fix is not a different
edge filter but a test: compare the observed difference between the states
with the difference an arbitrary split of the same samples produces.

### Added

- **Rewiring null** (`<comparison>/rewiring_vs_chance.csv`, and `q_rewiring`
  per metabolite in `rewired_metabolites.csv`). State labels are permuted
  across the pooled samples, both networks are rebuilt with the thresholds
  their observed counterparts used, and state-specific edges are counted. Each
  feature is centred within its state before permuting; without that the null
  is invalid, because mixing states that differ in abundance manufactures
  correlations inside every pseudo-group. Reports observed / null, an
  empirical p for the comparison, and an empirical p and BH q per feature.
- **Edge null** per state network (`details/network_*/null_model.csv`): each feature
  permuted independently within the state. Asks whether the network has more
  edges than features without joint structure would give.
- **Modularity against a configuration model** (same file): unweighted
  Louvain Q minus its mean over degree-preserving rewirings. Raw modularity
  largely restates density (Spearman between edge count and Q was −0.92
  across 30 real networks); delta Q does not.
- **Operating range.** Each comparison is classified as `within` (≥ 14
  samples in the smaller state), `reduced` (10–13) or `below` (< 10), printed
  by the input check and reported in `REPORT.md` and `report.json`. Runs below
  the range still execute, with a prominent warning.
- **Connectivity-based selection**, automatic. When every state of a comparison
  has at least `CONNECTIVITY.min_n_per_state` (14) samples, metabolites are also
  selected by change in their co-variation partners, whatever their abundance:
  the mean of each metabolite's 10 largest changes in Spearman rho, tested
  against permuted state labels with a pooled standardised null and
  Benjamini–Hochberg. Every selected metabolite carries `selection_route`
  (abundance / connectivity / both); the connectivity-only ones are what a
  differential-abundance analysis cannot report. Below 14 per state selection
  is by abundance only, with the reason in the log and the report. On simulated
  rewiring without any abundance change it selects 67% / 92% / 100% of the
  rewired metabolites at 14 / 20 / 25 samples per state (differential abundance:
  0%), with no false comparison-level call in 240 null and abundance-only
  replicates.
- **DGCA** (McKenzie et al., *BMC Syst Biol* 2016) as a second comparator in
  `dev/benchmarks/null_benchmark.py`, beside DiffCorr, with its R adapter and its
  rows in the committed result tables. It behaves as DiffCorr does: at 5-6 samples
  per state it calls a difference in every replicate, including when nothing
  differs. The DiffCorr and PhenoRewire rows are unchanged.
- **`dev/benchmarks/effect_size_benchmark.py`**, with its result tables: how large
  a reorganisation has to be to be detected at 5-12 samples per state, below the
  operating range. Pure rewiring with identical marginals is not detected there at
  all, since connectivity-based selection does not run below 14 per state.
- **Input check** (`checks.py`), at the start of every run and alone with
  `phenorewire --config ... --check`. Reads both tables the way the run will
  and reports features and samples, which samples matched, the groups and
  timepoints, and each comparison's sample size against the operating range.
  Stops with a message that says how to fix it when the data cannot be
  analysed: no matching samples, more than 1% negative intensities (autoscaled
  or centred data, which median normalisation and log2(x + 1) would silently
  corrupt), a missing timepoint, a group column that does not exist. Warns on
  features that hold text (binary variables such as `VRAI`/`FAUX`, previously
  dropped without notice) and on tables that look already log-transformed.
  Written to `run_info/data_check.txt`.
- **Readable terminal output.** Numbered steps with one or two lines of result
  each, warnings marked `!`, and a closing block with the result of each
  comparison and the files to open. Errors from the config or the data are
  printed as one message without a traceback (`--log-level DEBUG` shows it).
  The full log, with timestamps, is always written to `run_info/run.log`.
- `REPORT.md`: the result of each comparison in one sentence, warnings, the
  files to open, a table of the most rewired metabolites (name, rewired
  correlations, q, route, direction), one line per network, the methods and the
  reproducibility block. Replaces `final_report.md`.
- `rewired_metabolites.csv`: one row per selected metabolite with its name,
  selection route, direction, rewiring counts and per-metabolite test against
  chance, sorted by `q_rewiring` and then by `rewiring_score`, with that order
  carried in a `rank` column (v0.1.x sorted by `rewiring_score` alone).
  Replaces `priority_rewired_nodes.csv` and
  `rewiring_node_null.csv`, which had to be joined by hand. The same columns are
  node attributes of `rewiring_network.graphml`.
- `run_info/config_used.yaml`: the configuration as applied, every default
  filled in.
- **`NORM_METHOD`** accepts `pqn` (probabilistic quotient normalisation, for
  urine and other samples of variable dilution, and for sparse data), `tic`
  (total signal) and `none` (data already normalised outside PhenoRewire), in
  addition to `median`, which stays the default and is unchanged. Median
  normalisation counts undetected values, so a sample in which more than half
  of the features are absent has median 0 and was left unscaled; the input
  check now says so and points to `pqn`. Which method to choose, by sample type
  and number of features: `docs/prepare_your_data.md`. The method applied is
  printed by the input check and recorded in `REPORT.md`.
- Docs: `docs/prepare_your_data.md` (what the tables must contain; binary,
  scaled and log data; paired designs), `docs/reading_the_results.md` (every
  output file and column, for both modes; replaces `interpret_phenotype.md` and
  `interpret_temporal.md`, which now point to it), and `docs/configuration.md`
  (every configuration key, with a minimal config that the test suite validates
  against the bundled data).
- `dev/benchmarks/null_benchmark.py`: ground-truth benchmark of the released
  code, including the rewiring null, with DiffCorr as a comparator through R.

### Changed

- Selection permutations are drawn from one generator consumed feature by
  feature, the scheme of the code every published PhenoRewire analysis was run
  with (v0.1.x used one generator per feature). With the same seed, abundance-
  based selection reproduces those analyses exactly. Both schemes are valid
  Monte Carlo estimates of the same p-values; they differ only for features
  whose q-value sits at the threshold.
- Default `N_PERMUTATIONS_PHENO` and `TEMPORAL_CORRELATION.n_permutations`
  raised from 500 to 1000. At 500 the smallest attainable p (1/501) approaches
  the BH threshold for the leading features in a typical run. Configs that set
  the value explicitly are unaffected.
- The edge statistics (Spearman, BH, threshold mask) are factored out of
  `correlation_network()` into `pairwise_spearman()` and `edge_mask()`, so
  observed and null networks pass through one code path.
- README rewritten around what the tool measures and what it does not, with the
  benchmark behind each claim. Abundance-based selection is stated for what it
  is (a rank-based two-group test, equivalent to a permutation Mann–Whitney).
  Installation and the first run come first.
- **Output layout.** One folder per comparison, four files at its top, the
  rest under `details/` and `run_info/`, and no file written twice (a phenotype
  run wrote 79 files, 21 of them under a duplicated name; it now writes 30).
  All statistics are unchanged. Scripts that read 0.1.x outputs need the new
  paths:

  | 0.1.x | 0.2.0 |
  | --- | --- |
  | `report/final_report.md`, `triage/final_report.md` | `REPORT.md` |
  | `report/final_report.json` | `report.json` |
  | `triage/rewiring_null.csv` | `<comparison>/rewiring_vs_chance.csv` |
  | `triage/priority_rewired_nodes.csv`, `triage/rewiring_node_null.csv`, `rewiring_*/rewiring_node_summary.csv` | `<comparison>/rewired_metabolites.csv` |
  | `*_selection/selected_*_features_with_direction.csv`, `*_features_selected.csv` | `<comparison>/selected_metabolites.csv` |
  | `triage/rewiring_network_phenotype.graphml`, `..._temporal.graphml` | `<comparison>/rewiring_network.graphml` |
  | `*_selection/*_features_all.csv`, `connectivity_features_all.csv`, `*_selection_summary.json` | `<comparison>/details/selection/all_features_tested.csv`, `all_features_connectivity.csv`, `selection_summary.json` |
  | `network_<state>/nodes.csv`, `node_metrics.csv`, `triage/node_priority_ranking.csv`, `top_hubs.csv`, `top_brokers.csv`, `top_influencers.csv`, `priority_features.csv` | `<comparison>/details/network_<state>/nodes.csv` (all metrics and `triage_score`, best first) |
  | `network_<state>/community_summary.csv` | `<comparison>/details/network_<state>/communities.csv` |
  | `rewiring_*/edges_union_with_stats.csv`, `edges_shared.csv`, `edges_*_only.csv` | `<comparison>/details/rewiring/edges.csv` (column `state`) |
  | `rewiring_*/rewiring_summary.csv` | `<comparison>/details/rewiring/rewiring_summary.csv` |
  | `rewiring_*/metadata_used.csv` | `<comparison>/details/samples_used.csv` |
  | `metadata_aligned.csv`, `X_corr_matrix.csv`, `feature_annotation.csv` | `run_info/samples.csv`, `preprocessed_matrix.csv`, `feature_annotation.csv` |
  | `2dnetwork/`, `plots/` | `<comparison>/molecular_network/`, `<comparison>/plots/` |
  | `phenotype_direction.csv`, `group_counts_used.csv` | `report.json` and `run_info/data_check.txt` |

  `<comparison>` is `<ref>_vs_<case>` or `<group>_<T1>_vs_<T2>`. With
  `ANALYSIS_MODE: both`, `--run-2dnetwork` and `--plot-networks` now process
  both comparisons.
- `run._parse_time_numeric()` and `run._timepoint_label()` moved to
  `utils.parse_time_numeric()` and `utils.timepoint_label()`, unchanged, so the
  input check can read a design without importing the runner. Code that imported
  them from `phenorewire.run` must import them from `phenorewire.utils`.
- `triage.py` returns tables instead of writing them: `rank_network_nodes()`
  and `rewired_metabolites()` replace `export_network_triage()` and
  `export_rewiring_triage()`.
- A missing input file in the config is reported as a validation error naming
  the key (`DATA_MATRIX: file not found: ...`) instead of a `FileNotFoundError`
  traceback.

### Fixed

- The input check never reported samples that match no group in
  `GROUP_DEFINITION`: it tested the group labels for nullness, while
  `assign_groups()` labels an unmatched sample `"Unknown"`. QC injections, blanks
  and mistyped group values were dropped from a comparison without a word. They
  are now listed as a warning, by name.
- `compute_network_metrics()` dropped every node when given a feature table
  with rows but no columns, because pandas reports such a frame as empty.
  Unreachable from the CLI, where the table always carries columns.
- `tutorial_2dnetwork.ipynb` referred to an undefined `pr_path` and unpacked
  two of the three tables `build_annotation_expansion()` returns.

---

## [0.1.1] — 2026-09-02

Hardening and performance release. No change to the analysis method; results on
well-formed input are unchanged, except where noted under Fixed.

### Fixed

- **A requested analysis could be skipped in silence.** `ANALYSIS_MODE: phenotype`
  with three or more populated groups logged an INFO line, ran nothing, and exited 0
  with a report that claimed success and listed four output files that were never
  written. The mode is now resolved explicitly and raises with the group counts and
  the config keys to change. `PHENO_GROUP_REF`/`PHENO_GROUP_CASE`, when set, now
  restrict the run to that contrast instead of refusing because a third group exists —
  the config already validated those keys, but the runtime ignored them.
- **Spurious correlations from constant features.** Features with identical intensity
  in every sample (typically imputed at a fixed LOD value) were filtered only *after*
  median normalization. Column-wise scaling makes such a row proportional to the
  scale-factor vector, so it survived the variance check, correlated at rho = 1.0 with
  every other constant feature, and inherited whatever association the scale factors
  had with the phenotype. On a synthetic set with four imputed features, all 6 resulting
  edges were this artefact and all four features were selected as the *most*
  phenotype-associated in the dataset (p = 0.0005), ahead of every real one.

  Constant features are now identified on the raw intensities but dropped *after*
  normalization, so the scale factors are unchanged and every surviving feature keeps
  bit-identical values, r and p. Only features carrying no sample-to-sample information
  are removed; nothing else in the run moves.
- **Numeric metadata never matched the config.** `assign_groups()` compared numpy
  scalars from pandas against plain Python values, so `np.int64(1)` did not match the
  `1` written in the YAML. Affected metadata whose rows are not object-dtype.
- **Non-reproducible GraphML output.** Per-feature subnetwork exports iterated a Python
  set, so node and edge order followed per-process hash randomization and the same input
  produced byte-different files on consecutive runs.
- Empirical p-values depended on last-ulp float noise when permuted correlations tied
  exactly with the observed one — common with a discrete target. A tie tolerance makes
  them stable across BLAS builds and SciPy versions.
- `RANDOM_SEED: 0` was rejected by a shared ">= 1" validator.
- Config validation created `OUTDIR` as a side effect, so a config that failed on a
  later field still left an empty directory behind.
- `save_json()` wrote bare `NaN`/`Infinity` tokens, which strict JSON parsers reject,
  and could abort the whole run on an unencodable numpy scalar in the final step.
- Sub-hour timepoints collapsed onto the same directory name (`int(0.5)` and `int(0.8)`
  both gave `0h`). Fractional timepoints now render as `0p5h`.
- Duplicate or empty `FEATURE_ID_COL` values were accepted and produced self-loops and
  duplicated rows on every `.loc[]` lookup downstream; they are now rejected.
- `python -m phenorewire` silently lacked `--run-2dnetwork` and `--plot-networks`:
  `__main__.py` pointed at a second, less capable argparse entry point in `run.py`.
- The CLI's "install phenorewire[viz]" hint was unreachable, because `viz.py` imported
  matplotlib lazily inside the function; the user got a logged error instead.
- `correlation_network()` now warns when non-finite correlations discard edges. A single
  NaN removes a feature from the network entirely, previously without a word.

### Performance

- The Spearman permutation test is vectorized: ranking once up front turns each
  permutation into a matrix product instead of a fresh `scipy.stats.spearmanr` call.
  Measured at **272–338x** faster; the bundled demo run goes from 27s to under 6s.
  Verified bit-identical p-values against the previous implementation, with observed
  rho differing by at most 2.2e-16.

### Changed

- The phenotype and temporal selections share `selection.run_permutation_selection()`.
  They were ~90% identical and had drifted: only the phenotype path normalized the
  annotation `name` column, so temporal outputs could carry a missing or literal-"nan"
  name. Temporal outputs now get the same treatment.
- Removed `joblib` and `tqdm` from the runtime dependencies. `tqdm` was never imported;
  `joblib` is unnecessary after vectorization.
- `requirements.txt`, `requirements-dev.txt`, and `environment.yml` install the package
  rather than restating its dependencies, so the four lists cannot drift apart.
- Dropped the `fast_spearman_permutation_test` alias, `viz.plot_annotation_expansion_summary`,
  `triage._move_if_exists`, and a second argparse `main()` in `run.py` — all unused.
- k-hop neighbourhoods in the 2D annotation module are computed once instead of twice,
  and subnetwork filenames are sanitized against characters Windows rejects in paths.

### Added

- Ruff lint gate in CI, plus CI jobs for `python -m phenorewire` and for a
  runtime-dependencies-only install.
- Test coverage for `preprocessing.py`, `assign_groups()`, and run-mode resolution —
  previously at zero. The suite goes from 51 to 77 tests.

---

## [0.1.0] — 2026-06-09

### Summary
First public release. Intended for research use.

### Added

#### Package
- **Package renamed** `phenotype` → `phenorewire` for PyPI alignment
- `phenorewire/cli.py`: new CLI entry point with `--run-2dnetwork` and `--plot-networks` flags
- `phenorewire/__version__ = "0.1.0"`

#### 2D Network Module (`phenorewire/twodn/`)
- `chemical_integration.py`: `build_annotation_expansion()` — k-hop BFS in SpecReboot/GNPS
  molecular networks, filtered by `min_cosine`; produces TABLE 1 (summary) and TABLE 2 (expansion)
- `viz.py`: `plot_rewiring_network()` (PNG + optional HTML via pyvis) and
  `plot_annotation_expansion_summary()`
- New config parameters: `SPECREBOOT_NETWORK`, `PHENOREWIRE_NETWORK`, `MN_HOP_DEPTH` (2),
  `MN_MIN_COSINE` (0.70), `ANNOTATION_TOP_N` (20)
- Synthetic example GraphML files: `examples/data/example_specreboot_network.graphml` (30 nodes)
  and `examples/data/example_phenorewire_network.graphml` (22 nodes)

#### Parallelisation
- `spearman_permutation_test()` now uses `joblib.Parallel(n_jobs=-1, prefer="threads")` across
  features; per-feature seed = `(global_seed + feature_index) % 2**31` for full reproducibility
- Falls back to sequential loop if joblib is not installed
- `joblib>=1.2` added to core dependencies

#### Hardening
- Zero-variance features removed in preprocessing with a logged warning
- `TIMEPOINT_T1 < TIMEPOINT_T2` validated at config load time (raises `ValueError` if equal or reversed)
- Temporal analysis logs when other timepoints are present but silently excluded
- `pyproject.toml`: added `[all]` extra, `pyvis` to `[viz]`, `nbconvert` to `[dev]`

#### Tests
- 9 new unit tests in `tests/unit/test_2dnetwork.py`:
  1-hop correct, 2-hop correct, min_cosine filter, singleton handling, feature-not-in-MN,
  annotation_confidence, output files exist, all required summary columns, all required expansion columns
- Total: 51 tests (42 previous + 9 new), all green in < 30s

#### Repository structure
- Tests reorganised: `tests/unit/`, `tests/integration/`, `tests/simulate/`, `tests/benchmark/`
- `example/` moved to `examples/` with `examples/data/` subdirectory
- `examples/generate_example_graphml.py` generates synthetic GraphML files

#### Documentation
- Logo + PyPI/license/python badges added to README top
- `--run-2dnetwork` and `--plot-networks` section added to README
- New 2D network config parameters documented in README parameter reference
- Code structure section updated for `phenorewire/` and `phenorewire/twodn/`
- `python -m phenotype` corrected to `python -m phenorewire`
- `docs/interpretation_guide.md`: plain-language guide for biologists
- `CHANGELOG.md` (this file)
- FAQ: 2D network question added; test count updated to 51

### Changed

- Weight sensitivity simulation (`tests/simulate/weight_sensitivity.py`) redesigned:
  - Hub GT features (F000–F003) and hub BG decoys (F008–F012) share the same correlation block
  - Selection GT features (F004–F007) have strong mean shift but sparse connections
  - AUROC range ≥ 0.10 across weight pairs confirmed (topology-heavy best for hub patterns)
  - Previous design had SNR too high (AUROC=1.0 for all weight pairs — uninformative)

### Fixed
- `tests/benchmark/rewiring_benchmark.py`: `int()` cast prevents `ValueError: Unknown format code 'd' for float`

---

## [0.x] — Pre-release development

Internal development versions not tracked in this changelog.
See `git log` for detailed history.
