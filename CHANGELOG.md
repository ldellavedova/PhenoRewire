# Changelog

All notable changes to PhenoRewire are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
