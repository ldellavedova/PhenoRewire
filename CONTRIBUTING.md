# Contributing to PhenoRewire

Thank you for your interest in contributing. Contributions of any kind are welcome: bug reports, feature requests, documentation improvements, and code changes.

---

## Reporting issues

Use the [GitHub issue tracker](https://github.com/ldellavedova/PhenoRewire/issues) to report bugs, ask questions, or suggest features.

When reporting a bug, please include:
- your operating system and Python version
- the PhenoRewire version or commit hash (`pip show phenorewire`)
- the command you ran and the config file used (remove any private data)
- the relevant log output — re-run with `--log-level DEBUG` for more detail

---

## Development setup

```bash
git clone https://github.com/ldellavedova/PhenoRewire
cd PhenoRewire
conda env create -f environment.yml
conda activate phenorewire
pip install -e ".[dev]"
pytest dev/tests/
```

---

## Contributing code

1. Fork the repository and create a branch for your changes.
2. Install in development mode: `pip install -e ".[dev]"`
3. Make your changes and add tests in `dev/tests/` where appropriate.
4. Run the full test suite before submitting: `pytest dev/tests/`
5. Open a pull request with a short description of what you changed and why.

**Guidelines:**
- Keep changes focused: one fix or feature per pull request.
- Do not change the public CLI interface or config key names without first opening an issue for discussion.
- All new features must pass the existing test suite without regressions.
- Use clear commit messages.

---

## Running the example to verify changes

```bash
phenorewire --config config_phenotype.yaml
phenorewire --config config_phenotype.yaml --run-2dnetwork
phenorewire --config config_temporal.yaml
```

These commands use the bundled synthetic data in `data/` and should complete without errors.
