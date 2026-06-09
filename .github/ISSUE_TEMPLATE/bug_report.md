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

What happened instead? If there was an error, paste the full traceback below.
Re-run with `--log-level DEBUG` for more detail.

```
paste error or log output here
```

## Input data (if shareable)

Describe the size of your dataset (number of samples, features). If you can share a minimal reproducible example using the bundled demo data, please include the command.
