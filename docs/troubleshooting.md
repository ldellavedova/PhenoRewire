# PhenoRewire — Troubleshooting

---

## `phenorewire` command not found after installation

Activate your environment before running:
```bash
conda activate phenorewire   # or: source .venv/bin/activate
pip install -e .
```

---

## Python version mismatch

Use the Conda setup from `environment.yml`, which pins Python 3.11. Testing is run on Python 3.9–3.12; if you encounter version-specific issues, report them as a bug.

---

## Sample IDs do not match metadata

**Symptom:** metadata alignment fails, samples disappear, or group sizes are unexpectedly small.

**Fix:** ensure `META_SAMPLE_COL` values match the intensity-column names after stripping any prefix captured by `INTENSITY_REGEX`.

---

## Too few selected features

**Symptom:** the run stops at the feature selection step.

**Fix:** check `FDR_ALPHA`, `ADAPTIVE_SELECTION_THRESHOLDS`, `SELECTION_MIN_FEATURES`, and `MIN_FEATURES_HARD_STOP`. For small datasets, setting `MIN_FEATURES_HARD_STOP: 2` allows the run to proceed.

---

## Sparse networks

**Symptom:** the run stops after network construction.

**Fix:** check `CORR_ABS_THRESHOLD`, `CORR_FDR_ALPHA`, and the `NETWORK_MIN_*_HARD_STOP` parameters. For small demo datasets, `NETWORK_MIN_EDGES_HARD_STOP: 1` is appropriate.

---

## `--plot-networks` or the standalone figure script fails with ImportError

Install the plotting extras:
```bash
pip install "phenorewire[viz]"
```

---

## Molecular network node IDs do not match PhenoRewire feature IDs

**Symptom:** 2D outputs contain zero neighbours or many `in_specreboot = False` rows.

**Fix:** the molecular-network GraphML must use the same feature IDs as the rewiring GraphML for the features you want to expand.

---

## Missing cosine or similarity edge attributes

**Symptom:** 2D propagation expands through all edges regardless of confidence.

**Fix:** ensure the molecular-network GraphML stores a numeric similarity attribute such as `cosine`.

---

## Verify the installation

```bash
phenorewire --help
python -m pytest -q dev/tests/
```
