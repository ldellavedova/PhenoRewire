# PhenoRewire — Tutorials

This folder contains step-by-step tutorials for learning PhenoRewire.

---

## Tutorials

### [`tutorial.ipynb`](tutorial.ipynb) — Main phenotype pipeline walkthrough

Runs a complete PhenoRewire analysis on a synthetic dataset with known ground truth. Covers:
- exploring the input data
- configuring and running the pipeline via the Python API
- interpreting rewiring scores, sign-switch edges, and triage rankings
- comparing results to the embedded ground truth

**Best starting point** for new users who want to understand how the pipeline works.

### [`tutorial_2dnetwork.ipynb`](tutorial_2dnetwork.ipynb) — 2D molecular network annotation

Demonstrates the `--run-2dnetwork` feature using the bundled demo molecular network. Covers:
- running the annotation expansion via the Python API
- interpreting `annotation_summary.csv` and `annotation_expansion.csv`
- configuring hop depth and cosine threshold

**Prerequisite:** run `tutorial.ipynb` or the CLI quick-start first.

### [`walkthrough.md`](walkthrough.md) — Output interpretation walkthrough

A written step-by-step guide through the example outputs from the phenotype quick-start. Covers:
- reading `final_report.md`
- interpreting `priority_rewired_nodes.csv`
- understanding the network collapse in the example dataset
- visualising in Cytoscape

---

## Running the tutorials

Install PhenoRewire, then open a notebook:

```bash
pip install -e ".[dev]"
jupyter notebook docs/tutorials/tutorial.ipynb
```

The notebooks use the bundled data in `data/` and do not require additional downloads.
