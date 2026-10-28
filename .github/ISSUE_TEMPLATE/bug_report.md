---
name: Bug report
about: Report unexpected errors or incorrect output
title: "[BUG] "
labels: bug
assignees: ''
---

## Environment

- **OS:** (e.g., Ubuntu 22.04, macOS 14, Windows 11)
- **Python version:** (e.g., 3.11.4)
- **PhenoRewire version or commit:** (run `pip show phenorewire` or `git log -1 --format=%h`)
- **Installation method:** conda / venv / other

## Command run

```bash
phenorewire --config your_config.yaml [other flags]
```

## Config file

Paste the relevant parts of your config, replacing any private paths or group names:

```yaml
# paste relevant config here
```

## Expected behaviour

What did you expect to happen?

## Actual behaviour

What happened instead? Paste the message below, and attach `OUTDIR/run_info/run.log` and
`OUTDIR/run_info/data_check.txt` from the run. Re-run with `--log-level DEBUG` to see the full traceback.

```
paste error or log output here
```

## Input data (if shareable)

Describe the size of your dataset (number of samples, features, groups, timepoints): `data_check.txt` has all of it. Please do not attach study data. If the problem also happens with the bundled examples (`config_phenotype.yaml`, `config_temporal.yaml`), include the command.
