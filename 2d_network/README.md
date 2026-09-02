# 2D Network Module

This directory contains the standalone figure-generation workflow for the optional 2D molecular network view.

There are two related but different 2D workflows in this repository:

1. `phenorewire --config ... --run-2dnetwork`
   - package-integrated
   - produces CSV tables for neighbourhood-based annotation propagation over a molecular network
   - does not require users to manually inspect GraphML files

2. `python 2d_network/multilayer_2d_notebook.py --config ...`
   - standalone
   - produces the publication-style two-layer figure
   - is useful after you already have the PhenoRewire rewiring GraphML and a molecular network GraphML

> **The two workflows take different config files, and they are not interchangeable.**
>
> | File | Used by | Describes |
> |------|---------|-----------|
> | `config_phenotype.yaml`, `config_temporal.yaml` (repo root) | `phenorewire --config ...` | the analysis: input tables, groups, thresholds |
> | `2d_network/2d_figure_phenotype.yaml`, `2d_network/2d_figure_temporal.yaml` | `multilayer_2d_notebook.py --config ...` | the figure: layout, colours, labels, output format |
>
> Passing one where the other is expected will fail with a validation error, not
> produce a wrong result.

## Package-integrated annotation workflow

Use the main CLI when you want CSV outputs such as:

- `annotation_summary.csv`
- `annotation_expansion.csv`
- `annotation_first_neighbours.csv`

Example:

```bash
phenorewire --config config_phenotype.yaml --run-2dnetwork
```

Required config keys in the main PhenoRewire config:

- `SPECREBOOT_NETWORK`
- `PHENOREWIRE_NETWORK` (optional; auto-detected from the current run if omitted)
- `MN_HOP_DEPTH`
- `MN_MIN_COSINE`
- `ANNOTATION_TOP_N`

## Standalone figure workflow

Install plotting extras:

```bash
pip install "phenorewire[viz]"
```

Then run from the repository root:

```bash
python 2d_network/multilayer_2d_notebook.py --config 2d_network/2d_figure_phenotype.yaml
```

The standalone script expects:

- a PhenoRewire rewiring GraphML
- a molecular/chemical similarity GraphML
- an annotation table for labelling the lower network layer

See `README.md` for the main user-facing installation and workflow instructions.
