"""
Generate the two bundled example datasets in data/.  Both are synthetic.

Run from the repository root:
    python dev/benchmarks/fixtures/generate_example_data.py

Phenotype example (feature_table.csv, metadata.csv, config_phenotype.yaml):
60 metabolites x 40 samples (20 reference, 20 case).  Twenty per state puts
the example inside PhenoRewire's operating range, so a run shows every part of
the workflow: selection by abundance and by co-variation partners, the
networks, the rewiring classification and the null that tests it.

Each block of metabolites is one case from the benchmark
(dev/benchmarks/null_benchmark.py), so the example has a known answer:

  F001-F010  abundance + rewiring   correlated in reference, independent in
                                    case, and higher in case.  Selected by
                                    abundance; their edges are reference-only.
  F011-F020  abundance only         correlated in both states, lower in case.
                                    Selected by abundance; their edges are
                                    shared, so they are not rewired.
  F021-F030  rewiring only          independent in reference, correlated in
                                    case, with identical means and spreads.
                                    Invisible to differential abundance;
                                    selected by co-variation partners
                                    (selection_route = connectivity); their
                                    edges are case-only.
  F031-F060  background             nothing; should not be selected.

One block goes up and one goes down, and half the table is background, so
sample-wise median normalisation leaves the shifts in place, as it does with
real data.  The names are illustrative labels for a demo, not a biological
model; F001-F012 match the nodes of data/example_molecular_network.graphml.

Temporal example (temporal_feature_table.csv, temporal_metadata.csv,
config_temporal.yaml): the same 60 metabolites in one group of 20 subjects,
each sampled at 0, 4 and 10 h (60 samples).  The networks compare 0 h with
10 h; selection uses all three timepoints, so 4 h sits halfway in every
change.  The blocks mirror the phenotype example:

  F001-F010  abundance + rewiring   correlated at 0 h, independent at 10 h,
                                    rising over time.
  F011-F020  abundance only         correlated throughout, falling over time.
  F021-F030  rewiring only          independent at 0 h, correlated at 10 h,
                                    same means and spreads at every timepoint.
  F031-F060  background             nothing.

The subject column records the repeated sampling, as a real design would; the
values are drawn independently per sample, and PhenoRewire does not use the
pairing (see docs/prepare_your_data.md).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
TEMPORAL_SEED = 7
N_PER_STATE = 20
TIMEPOINTS = (0, 4, 10)
N_FEATURES = 60
BLOCK = 10
OUT = Path(__file__).resolve().parents[3] / "data"

ANNOTATIONS = [
    ("F001", 133.014, 0.98, "lactate"),
    ("F002", 175.119, 2.31, "citrulline"),
    ("F003", 146.058, 1.87, "glutamate"),
    ("F004", 132.054, 3.12, "asparagine"),
    ("F005", 204.090, 4.55, "adenosine"),
    ("F006", 267.097, 5.21, "AMP"),
    ("F007", 348.069, 4.83, "ADP"),
    ("F008", 428.037, 4.91, "ATP"),
    ("F009", 117.043, 1.44, "threonine"),
    ("F010", 131.058, 1.69, "methionine"),
    ("F011", 147.053, 1.55, "glutamine"),
    ("F012", 132.102, 2.04, "leucine"),
    ("F013", 168.065, 3.77, "phenylalanine"),
    ("F014", 182.081, 4.22, "tyrosine"),
    ("F015", 204.090, 5.01, "tryptophan"),
    ("F016", 191.020, 1.11, "citrate"),
    ("F017", 173.009, 0.95, "isocitrate"),
    ("F018", 90.032, 0.72, "pyruvate"),
    ("F019", 116.047, 1.33, "succinate"),
    ("F020", 245.077, 6.14, None),
    ("F021", 309.045, 7.88, None),
    ("F022", 181.061, 8.33, None),
] + [(f"F{i:03d}", round(100.0 + 7.3 * i, 3), round(0.5 + 0.13 * i, 2), None) for i in range(23, 61)]


def block(rng, n: int, k: int, rho: float) -> np.ndarray:
    """k features x n samples, pairwise correlation rho."""
    if rho <= 0:
        return rng.standard_normal((k, n))
    return np.sqrt(rho) * rng.standard_normal((1, n)) + np.sqrt(1 - rho) * rng.standard_normal((k, n))


def standardise(m: np.ndarray) -> np.ndarray:
    return (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, ddof=1, keepdims=True)


def write_table(X: np.ndarray, samples: list[str], path: Path) -> None:
    """log2-scale values -> intensities, with the annotation columns."""
    baseline = np.linspace(6.0, 9.0, N_FEATURES)[:, None]   # log2 scale, per metabolite
    intensities = np.round(2.0 ** (X * 0.6 + baseline), 4)
    rows = []
    for (fid, mz, rt, name), values in zip(ANNOTATIONS, intensities):
        row = {"feature_id": fid, "mz": mz, "rt": rt, "name": name or ""}
        row.update(dict(zip(samples, values)))
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def phenotype_example() -> None:
    rng = np.random.default_rng(SEED)
    n = N_PER_STATE
    ref = rng.standard_normal((N_FEATURES, n))
    case = rng.standard_normal((N_FEATURES, n))
    d, a, b = slice(0, BLOCK), slice(BLOCK, 2 * BLOCK), slice(2 * BLOCK, 3 * BLOCK)

    # F001-F010: correlated in reference, independent and higher in case.
    ref[d] = block(rng, n, BLOCK, 0.85)
    case[d] = block(rng, n, BLOCK, 0.0) + 1.5
    # F011-F020: correlated in both, lower in case.
    ref[a] = block(rng, n, BLOCK, 0.8)
    case[a] = block(rng, n, BLOCK, 0.8) - 1.5
    # F021-F030: independent in reference, correlated in case, same marginals.
    ref[b] = standardise(block(rng, n, BLOCK, 0.0))
    case[b] = standardise(block(rng, n, BLOCK, 0.85))

    samples = [f"Sample_{i:02d}" for i in range(1, 2 * n + 1)]
    OUT.mkdir(parents=True, exist_ok=True)
    write_table(np.hstack([ref, case]), samples, OUT / "feature_table.csv")
    pd.DataFrame({
        "sample_id": samples,
        "group_label": ["reference"] * n + ["case"] * n,
        # Half of each group at timepoint 1 and half at 2, so config_temporal.yaml
        # has a structure to run on; the phenotype example ignores it.
        "timepoint": ([1] * (n // 2) + [2] * (n - n // 2)) * 2,
    }).to_csv(OUT / "metadata.csv", index=False)
    print(f"wrote {OUT / 'feature_table.csv'} ({N_FEATURES} x {2 * n}) and metadata.csv")


def temporal_example() -> None:
    rng = np.random.default_rng(TEMPORAL_SEED)
    n = N_PER_STATE
    d, a, b = slice(0, BLOCK), slice(BLOCK, 2 * BLOCK), slice(2 * BLOCK, 3 * BLOCK)
    # Fraction of the way from 0 h to 10 h; 4 h sits halfway.
    progress = {0: 0.0, 4: 0.5, 10: 1.0}

    parts, samples, meta = [], [], []
    for t in TIMEPOINTS:
        f = progress[t]
        X = rng.standard_normal((N_FEATURES, n))
        # F001-F010: coupling fades and abundance rises over time.
        X[d] = block(rng, n, BLOCK, 0.85 * (1 - f)) + 1.5 * f
        # F011-F020: coupled throughout, abundance falls over time.
        X[a] = block(rng, n, BLOCK, 0.8) - 1.5 * f
        # F021-F030: coupling appears over time, marginals unchanged.
        X[b] = standardise(block(rng, n, BLOCK, 0.85 * f))
        parts.append(X)
        for s in range(1, n + 1):
            samples.append(f"S{s:02d}_{t}h")
            meta.append({"sample_id": f"S{s:02d}_{t}h", "subject": f"S{s:02d}",
                         "cohort": "exposed", "timepoint": t})

    write_table(np.hstack(parts), samples, OUT / "temporal_feature_table.csv")
    pd.DataFrame(meta).to_csv(OUT / "temporal_metadata.csv", index=False)
    print(f"wrote {OUT / 'temporal_feature_table.csv'} ({N_FEATURES} x {len(samples)}) and temporal_metadata.csv")


def main() -> None:
    phenotype_example()
    temporal_example()


if __name__ == "__main__":
    main()
