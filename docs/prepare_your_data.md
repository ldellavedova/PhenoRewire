# Preparing your data

PhenoRewire reads two comma-separated tables. This page says what they must contain, what to leave out, and how to check them before running.

```bash
phenorewire --config my_config.yaml --check
```

reads both tables exactly as a run would, and reports what it found and anything that has to be fixed, without running the analysis. Every run starts with the same check.

---

## The feature table

One row per metabolite feature; identifier columns, then one intensity column per sample.

| feature_id | mz | rt | name | Sample_01 | Sample_02 | … |
|---|---|---|---|---|---|---|
| F001 | 133.014 | 0.98 | lactate | 15230.4 | 18820.1 | … |
| F002 | 175.119 | 2.31 | citrulline | 820.7 | 0 | … |

- **`feature_id`** must be unique. A compound name usually is not.
- **`mz`, `rt`** are required; **`name`** is optional (set `NAME_COL`) and is used in every table and in the report.
- **Intensity columns** are picked by `INTENSITY_REGEX`. With MZmine-style headers (`Sample_01.mzML Peak area`) the suffix and the file extension are removed to obtain the sample name, which must match the metadata. See [`configuration.md`](configuration.md#intensity_regex).
- **Missing values**: leave the cell empty or write 0. Both mean "not detected". Features detected in fewer than `MIN_PRESENCE_GLOBAL` (default 25%) of samples are removed.

### Intensities: raw, non-negative, not scaled

PhenoRewire normalises each sample (`NORM_METHOD`, next section) and applies log2(x + 1) itself. It needs the intensities before scaling:

| Your table holds | What to do |
|---|---|
| Peak areas or heights from MZmine, MS-DIAL, XCMS, … | Use them as they are, and choose `NORM_METHOD` below |
| Intensities already normalised to an external factor (creatinine, specific gravity, osmolality, cell count, internal standards, PQN) | Use them with `NORM_METHOD: none`, so they are not normalised a second time |
| **Log-transformed** intensities | Set `LOG_TRANSFORM: false`, or provide the values before the log transform. The check warns when every value is below 50 |
| **Autoscaled, Pareto-scaled or mean-centred** values (negative numbers) | **Not usable.** Values ≤ 0 count as "not detected", and normalising centred data is meaningless. Provide the values before scaling. The check stops the run when more than 1% of values are negative |
| Imputed values (e.g. half the detection limit) | Usable; features that are constant in every sample are removed automatically |

### Choosing the normalisation (`NORM_METHOD`)

Normalisation removes differences between samples that have nothing to do with biology: how much material was extracted or injected, or how dilute a urine sample was. For a co-variation network this matters more than for a differential-abundance test. A dilution factor left in the data is shared by every metabolite of a sample, so every pair of metabolites correlates through it and the network describes dilution instead of metabolism. A normalisation that does not fit the data can create shared variation in the same way.

| Your data | `NORM_METHOD` | Why |
|---|---|---|
| Untargeted peak areas from samples of comparable concentration: plasma, serum, cells, tissue extracts, culture supernatants | `median` (default) or `pqn` | Corrects the amount extracted or injected. Both assume that most metabolites do not change between samples |
| Urine, faeces, saliva, sweat, or any sample whose dilution varies | `pqn` | The standard correction for dilution (Dieterle et al., *Anal Chem* 2006), and robust when many metabolites change |
| Sparse data: many features undetected in many samples (exposome, microbial secondary metabolites) | `pqn` | `median` counts undetected values as 0: a sample in which more than half of the features are absent has median 0 and is not scaled at all. The input check warns when this happens. `pqn` uses only the detected values |
| Already normalised outside PhenoRewire (see the table above), or absolute concentrations from a calibrated targeted method | `none` | The values are already comparable between samples |
| Replicating an analysis that used total-signal normalisation | `tic` | Only then. If a few abundant metabolites change between the states, `tic` passes that change on to every other metabolite of the sample |

**How many features.** `median` and `pqn` estimate one factor per sample from all its features, and assume that most of them do not change. With an untargeted table of hundreds of features that is a safe assumption. With a small targeted panel (a few dozen compounds, many of which may change), the estimate is noisy or biased: prefer an external normaliser (specific gravity, creatinine, internal standards) and `NORM_METHOD: none`.

**How many samples.** The choice does not depend on the number of samples per state. `pqn` builds its reference profile from all samples together; it works from about ten samples in total.

**How `pqn` is computed here.** The reference is each feature's median across samples, built from the features detected in every sample when there are at least ten of them (otherwise from all features, over the samples where each was detected). A sample's factor is the median of its intensities divided by the reference, over its detected features.

If unsure: `median` for blood, cells and tissue; `pqn` for urine and for sparse data; `none` when you normalised before.

The analyses in the paper used `median` on culture supernatants. `pqn` and `tic` are offered for other matrices and were not used there.

### Binary variables

Leave out variables that are only detected / not detected (0/1, TRUE/FALSE, or `VRAI`/`FAUX` from a French-language Excel). Correlation networks are built from continuous intensities: a Spearman correlation between two binary variables is dominated by ties and is not comparable to one between intensities. Text values cannot be read as intensities; the check names those features and they are removed. If detection itself is of interest, analyse it separately (for example with Fisher's or McNemar's test).

---

### Quality-control samples, blanks and pools

Pooled QC injections, procedural and solvent blanks, and standard mixes are not biological samples and must not enter a network. A pooled QC is, by construction, the average of the study: a run of QC injections carries almost no biological variation, so including them compresses the correlations of every metabolite towards the pool. Blanks add the opposite problem, a column of near-zero intensities that makes unrelated features look as though they co-vary.

**Remove their intensity columns from the feature table** before running. Giving them no group in `GROUP_DEFINITION` keeps them out of the comparison — the input check lists every sample it drops that way — but not out of the two steps that run before groups are assigned: the presence filter, which counts in how many of *all* samples a feature was detected, and the normalisation reference, which is taken over all columns. A handful of blanks can therefore push a real feature below `MIN_PRESENCE_GLOBAL`. Leaving them ungrouped does not change which correlations are found (a shared scale factor cannot), but it does change which features are tested.

Keep them for the work that belongs before PhenoRewire: blank subtraction or a blank-ratio filter, and a QC-based filter on the features whose relative standard deviation across the QCs is too high. PhenoRewire expects a table on which that filtering has already been done.

### Two ionisation modes

Run positive and negative mode separately, with their own config and `OUTDIR`, and compare the two result sets afterwards. A feature table is one measurement space: concatenating the modes puts the same metabolite in twice, under two identifiers, and the pair then correlates almost perfectly and enters every network as a spurious edge. Agreement between the modes is a useful check on a result, not something to average over.

---

## The metadata table

One row per sample.

| sample_id | condition | timepoint | subject |
|---|---|---|---|
| Sample_01 | control | 0 | P01 |
| Sample_02 | treated | 0 | P02 |

- **`META_SAMPLE_COL`** names the column that matches the sample names of the feature table. Every sample of the feature table needs a row; extra rows are ignored.
- **Groups** are defined in the config from any metadata columns:

  ```yaml
  GROUP_DEFINITION:
    Control: {condition: "control"}
    Treated: {condition: "treated", timepoint: 0}   # several conditions = all must match
  ```

  Values must match exactly (`control` is not `Control`). Samples that match no group are not used.
- **Time** (temporal mode) is a numeric column: `3`, `3h` and `3 hours` all read as 3.

---

## Choosing the comparison

| Design | Mode | Config |
|---|---|---|
| Two groups (case vs control) | `phenotype` | two entries in `GROUP_DEFINITION`, `PHENO_GROUP_REF`, `PHENO_GROUP_CASE` |
| Two or more timepoints within one group | `temporal` | `TEMPORAL_CORRELATION` (enabled, group, time column), `TIMEPOINT_T1`, `TIMEPOINT_T2` |
| Both, in one run | `both` | both of the above |

In temporal mode, metabolites are selected by their association with time across **all** timepoints of the group; the two networks are then built at `TIMEPOINT_T1` and `TIMEPOINT_T2`. To compare other pairs of timepoints, run again with a different pair and a different `OUTDIR`.

### Sample size

Each network is estimated from the samples of one state. The input check reports the number per state and whether it is within the operating range:

| Samples in the smaller state | Status |
|---|---|
| ≥ 14 | within the operating range; metabolites are also selected by change in co-variation partners |
| 10–13 | reduced: selection by abundance only, low recall |
| < 10 | below: treat any result as preliminary |

### Paired and repeated-measures designs

If the same subjects are sampled at both timepoints (or in both conditions), PhenoRewire does not use the pairing. Each network is built from one sample per subject, so the networks themselves are unaffected. The permutation tests treat the samples as unpaired; for positively correlated repeated measures this makes them conservative (valid, with somewhat less power). Keep a subject column in the metadata anyway: it costs nothing and keeps the design explicit.

### Pooled conditions

If a state pools several distinct conditions (several strains, sites or batches), metabolites that differ between those conditions correlate across the pooled samples, and that correlation is not co-regulation. Compare conditions separately, or correct for batch before running.
