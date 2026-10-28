# PhenoRewire: example walkthrough

A guided reading of one run on the bundled example data, in the order the
outputs should be read. Every number below comes from that run.

## The example data

| | |
| --- | --- |
| Metabolites | 60 |
| Samples | 40: 20 reference, 20 case |
| Comparison | phenotype mode, reference vs case |

The data are synthetic (`dev/benchmarks/fixtures/generate_example_data.py`)
and built so that the right answer is known. Each block of ten metabolites is
one of the situations the benchmark tests:

| Metabolites | What happens between reference and case | What a correct analysis reports |
| --- | --- | --- |
| F001–F010 | correlated in reference, independent in case, and higher in case | changed in abundance **and** rewired |
| F011–F020 | correlated in both states, lower in case | changed in abundance, **not** rewired |
| F021–F030 | independent in reference, correlated in case; level and spread identical | **rewired without any change in abundance**: invisible to differential abundance |
| F031–F060 | background | nothing |

The names on F001–F022 (lactate, citrulline, ATP...) are labels for the demo,
not a biological model.

## Step 1: check, then run

From the repository root:

```bash
phenorewire --config config_phenotype.yaml --check
phenorewire --config config_phenotype.yaml
```

The check reports 60 features x 40 samples, all samples matched, the two
groups with 20 samples each, and that the comparison is within the operating
range. Twenty samples per state means PhenoRewire selects metabolites both by
abundance and by change in co-variation partners.

The run takes a few seconds and ends with:

```text
Done in 8 s.
  Reference vs Case: the networks differ more than chance (86 state-specific correlations vs
  35.7 expected, p = 0.005); 18 metabolite(s) rewire more than chance (q <= 0.05).

  Start here : results/phenotype_example/REPORT.md
  Key tables : results/phenotype_example/Reference_vs_Case/
```

Open `REPORT.md`: the same sentence comes first, then the files to open and the
most rewired metabolites. The steps below go through those files, which are
all in `results/phenotype_example/Reference_vs_Case/`.

## Step 2: is there rewiring at all?

Open `rewiring_vs_chance.csv` first. It answers whether the two states differ
in co-variation by more than an arbitrary split of the same 40 samples:

| observed state-specific edges | expected by chance | ratio | p |
| --- | --- | --- | --- |
| 86 | 35.7 | 2.4 | 0.005 |

Two and a half times what chance gives, at the smallest p 200 permutations
can return. Only now is it worth reading further. Had the ratio been near 1,
the rankings below would describe noise.

## Step 3: which metabolites rewired?

Open `rewired_metabolites.csv`. One row per selected metabolite, sorted by
`q_rewiring`: whether it carries more state-specific edges than chance.

| rank | metabolites | `q_rewiring` | `selection_route` | what it is |
| --- | --- | --- | --- | --- |
| 1–10 | F001–F010 | 0.009 | both | lost their partners in case, and rose in abundance |
| 11–18 | F021–F026, F028, F029 | 0.009 | connectivity | gained partners in case; abundance unchanged |
| 19 | F013 | 0.09 | abundance | 7 state-specific edges against 3.8 expected: not significant |
| 20– | the rest | ≥ 0.20 | | |

The 18 metabolites with `more_than_chance = True` are exactly the rewired ones
in the design, and none of F011–F020 is among them. Their block has 26 shared
edges but also 15 state-specific ones that crossed the edge threshold in one
state and not the other by chance; the most affected, F013, stays above
q = 0.05. This is exactly the error a raw count of state-specific edges makes,
and what the null prevents. None of the background metabolites selected by
chance is significant either (best: F052, q = 0.20).

Within the significant metabolites, `rewiring_score` (state-specific plus
sign-switching edges) and `rewiring_fraction` say how much of each one's
neighbourhood changed.

## Step 4: which metabolites were selected, and how?

Open `selected_metabolites.csv`. 32 metabolites were selected; the column
`selection_route` says why:

| route | metabolites | meaning |
| --- | --- | --- |
| `both` | F001–F010 | changed in abundance and in co-variation partners |
| `abundance` | F011–F020, F038 | changed in abundance |
| `connectivity` | F021–F026, F028, F029, F050, F052, F058 | changed partners; abundance unchanged |

The `connectivity` route is what a differential-abundance analysis cannot give
you. Compare F021 with F001:

| metabolite | correlation with state | abundance q | log2 fold change | route | connectivity q |
| --- | --- | --- | --- | --- | --- |
| F001 | +0.68 | 0.003 | +1.12 | both | 0.012 |
| F021 | −0.07 | 0.91 | −0.06 | connectivity | 0.020 |

F021 has no abundance signal at all (q = 0.91), yet its partners changed.

Four background metabolites (F038, F050, F052, F058) passed selection by
chance, as expected at a 5% false-discovery rate; step 3 showed that none of
them rewires. The statistics of all 60 metabolites, selected or not, are in
`details/selection/`.

## Step 5: the two networks

`REPORT.md` has one line per network; the files are in
`details/network_Reference/` and `details/network_Case/`.

| network | nodes | edges | modularity | edges vs chance | modularity beyond degree sequence (ΔQ) |
| --- | --- | --- | --- | --- | --- |
| Reference | 32 | 81 | 0.50 | p = 0.005 | +0.33 |
| Case | 32 | 57 | 0.50 | p = 0.005 | +0.30 |

Both networks hold far more edges than metabolites without joint structure
would give, and both are modular beyond what their degree sequence imposes
(`null_model.csv`). Where the edges sit (`details/rewiring/edges.csv`, column
`state`):

| pair of blocks | reference only | case only | shared |
| --- | --- | --- | --- |
| F001–F010 among themselves | 43 | 0 | 0 |
| F011–F020 among themselves | 12 | 3 | 26 |
| F021–F030 among themselves | 0 | 27 | 0 |

The first block lost its co-variation in the case state, the third gained it,
and the second kept it.

## Step 6: visualise (optional)

Open `rewiring_network.graphml` in Cytoscape:

1. File → Import → Network from File.
2. Colour edges by `state` (`Reference_only`, `Case_only`, `shared`).
3. Size nodes by `rewiring_score`, and highlight nodes with `more_than_chance`.

## Step 7: annotation from a molecular network (optional)

```bash
phenorewire --config config_phenotype.yaml --run-2dnetwork
```

places the top rewired metabolites in the demo molecular network
(`data/example_molecular_network.graphml`) and suggests annotations from their
spectral neighbours; see `Reference_vs_Case/molecular_network/annotation_summary.csv`
and [`interpret_2d_network.md`](../interpret_2d_network.md).

## The temporal example

```bash
phenorewire --config config_temporal.yaml
```

runs the same analysis on a synthetic time course: one group of 20 subjects,
each sampled at 0, 4 and 10 h. Metabolites are selected by their association
with time across all three timepoints; the two networks compare 0 h with 10 h.
The blocks mirror the phenotype example (F001–F010 rise and lose their
partners, F011–F020 fall and keep them, F021–F030 gain partners without
changing in abundance), and the outputs are in
`results/temporal_example/Exposed_0h_vs_10h/`:

- `rewiring_vs_chance.csv`: 96 state-specific edges against 33.5 expected,
  p = 0.005.
- `rewired_metabolites.csv`: 17 metabolites rewire more than chance, nine of
  F001–F010 and eight of F021–F030, the latter selected only by connectivity.
  None of F011–F020 does, although their block has 30 edges present at one
  timepoint only: sampling noise, which the null recognises.

## Output files

```text
results/phenotype_example/
├── REPORT.md                        start here
├── Reference_vs_Case/
│   ├── rewiring_vs_chance.csv       step 2: is there rewiring at all?
│   ├── rewired_metabolites.csv      step 3: which metabolites rewired?
│   ├── selected_metabolites.csv     step 4: what was selected, and how
│   ├── rewiring_network.graphml     step 6: Cytoscape
│   └── details/                     selection statistics, each network, every edge
└── run_info/                        log, input check, configuration used
```

Every column is explained in [`reading_the_results.md`](../reading_the_results.md).
