# Changelog

All notable changes to PhenoRewire are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
