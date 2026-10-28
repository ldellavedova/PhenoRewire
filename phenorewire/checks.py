# ── PhenoRewire · Input check ────────────────────────────────────────────────
# Reads the two input tables the way the pipeline will, before any statistics,
# and says in plain words what it found: how many features and samples, which
# samples matched, the groups and timepoints, the comparisons that will run and
# whether their sample sizes are in the operating range.  Problems that would
# make the results meaningless (no matching samples, negative intensities, a
# missing timepoint) stop the run here with a message that says how to fix them.
#
# Runs at the start of every run, and alone with `phenorewire --config ... --check`.
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .preprocessing import intensity_columns, matrix_sample_ids, norm_sample_id
from .utils import NORM_LABELS, assign_groups, parse_time_numeric, timepoint_label

# Operating range for network inference, in samples per state.  Measured on
# simulated ground truth (dev/benchmarks/null_benchmark.py): below
# OPERATING_RANGE_MIN_N no method tested recovers rewiring reliably, and below
# OPERATING_RANGE_RECOMMENDED_N recall is low.  Runs outside the range still
# execute, because a result may be needed for comparison, but every output
# that depends on the networks carries the status.
OPERATING_RANGE_MIN_N = 10
OPERATING_RANGE_RECOMMENDED_N = 14

# Above this share of negative intensities the table cannot be raw intensities.
NEGATIVE_SHARE_ERROR = 0.01
# Raw LC-MS intensities run to thousands and more; a table whose largest value
# is below this has almost certainly been log-transformed already.
LOG_SCALE_MAX = 50.0


class InputDataError(ValueError):
    """The input tables cannot be analysed as they are."""


def operating_range(n_a: int, n_b: int) -> dict:
    """Where a two-state comparison sits relative to the measured operating range."""
    n_min = int(min(n_a, n_b))
    if n_min >= OPERATING_RANGE_RECOMMENDED_N:
        status = "within"
    elif n_min >= OPERATING_RANGE_MIN_N:
        status = "reduced"
    else:
        status = "below"
    return {
        "n_per_state": [int(n_a), int(n_b)],
        "min_n_per_state": n_min,
        "recommended_min_n": OPERATING_RANGE_RECOMMENDED_N,
        "absolute_min_n": OPERATING_RANGE_MIN_N,
        "status": status,
    }


@dataclass
class DataCheck:
    data_matrix: str = ""
    metadata: str = ""
    lines: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    comparisons: list[dict] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def text(self) -> str:
        out = [
            "PhenoRewire input check",
            f"  feature table : {self.data_matrix}",
            f"  metadata      : {self.metadata}",
            "",
        ]
        out += [f"  {line}" for line in self.lines]
        if self.warnings:
            out += ["", "Warnings"]
            out += [f"  ! {w}" for w in self.warnings]
        if self.errors:
            out += ["", "Problems that must be fixed before running"]
            out += [f"  x {e}" for e in self.errors]
        out += ["", "Result: " + ("ready to run." if self.ok else f"{len(self.errors)} problem(s) to fix.")]
        return "\n".join(out)

    def error_message(self) -> str:
        n = len(self.errors)
        head = f"The input data cannot be analysed as they are ({n} problem{'s' if n != 1 else ''}):"
        return "\n".join([head] + [f"- {e}" for e in self.errors])


def _examples(values, k: int = 5) -> str:
    vals = [str(v) for v in values]
    shown = ", ".join(vals[:k])
    return shown + (f", ... ({len(vals)} in total)" if len(vals) > k else "")


def _pct(x: float) -> str:
    return f"{100 * x:.1f}%" if x < 0.995 else f"{100 * x:.0f}%"


def _check_feature_table(config, df: pd.DataFrame, chk: DataCheck) -> Optional[list[str]]:
    """Checks the feature table; returns the intensity columns, or None when unusable."""
    required = {
        "FEATURE_ID_COL": config.FEATURE_ID_COL,
        "MZ_COL": config.MZ_COL,
        "RT_COL": config.RT_COL,
    }
    for key, col in required.items():
        if col not in df.columns:
            chk.errors.append(
                f"The feature table has no column '{col}' ({key} in the config). "
                f"Its first columns are: {_examples(df.columns, 8)}."
            )
    if config.NAME_COL and config.NAME_COL not in df.columns:
        chk.warnings.append(
            f"NAME_COL = '{config.NAME_COL}' is not a column of the feature table; "
            "metabolites will be shown by their feature ID."
        )

    cols = intensity_columns(df, config.INTENSITY_REGEX)
    if not cols:
        chk.errors.append(
            f"No column of the feature table matches INTENSITY_REGEX = '{config.INTENSITY_REGEX}'. "
            f"Its first columns are: {_examples(df.columns, 8)}."
        )
        return None
    chk.lines.append(f"{len(df)} features x {len(cols)} samples in the feature table")

    if config.FEATURE_ID_COL in df.columns:
        ids = df[config.FEATURE_ID_COL].astype(str).str.strip()
        dup = ids[ids.duplicated(keep=False)].unique().tolist()
        if dup:
            chk.errors.append(
                f"{len(dup)} feature ID(s) appear more than once in '{config.FEATURE_ID_COL}' "
                f"(e.g. {_examples(dup)}). Every feature needs a unique ID."
            )
        empty = int(ids.isin({"", "nan", "None"}).sum())
        if empty:
            chk.errors.append(f"{empty} feature(s) have no ID in '{config.FEATURE_ID_COL}'.")

    raw = df[cols]
    num = raw.apply(pd.to_numeric, errors="coerce")
    text_cells = raw.notna() & num.isna()
    all_text = pd.Series(False, index=df.index)
    if text_cells.values.any():
        values = pd.Series(raw.values[text_cells.values]).astype(str).str.strip()
        common = values.value_counts().index.tolist()
        all_text = text_cells.sum(axis=1).eq(raw.notna().sum(axis=1)) & text_cells.any(axis=1)
        some_text = text_cells.any(axis=1) & ~all_text
        if all_text.any():
            chk.warnings.append(
                f"{int(all_text.sum())} feature(s) hold text instead of intensities "
                f"(values such as {_examples([repr(v) for v in common], 3)}) and will be removed. "
                "If they are binary (detected / not detected) variables, remove them from the table: "
                "correlation networks need continuous intensities."
            )
        if some_text.any():
            chk.warnings.append(
                f"{int(some_text.sum())} feature(s) contain some text values "
                f"(such as {_examples([repr(v) for v in common], 3)}); those cells are treated as not detected."
            )

    vals = num.values.astype(float)
    finite = np.isfinite(vals)
    n_finite = int(finite.sum())
    if n_finite == 0:
        chk.errors.append("The intensity columns contain no numbers.")
        return cols

    neg = finite & (vals < 0)
    neg_share = neg.sum() / n_finite
    if neg_share > NEGATIVE_SHARE_ERROR:
        n_rows = int(neg.any(axis=1).sum())
        chk.errors.append(
            f"{_pct(neg_share)} of the intensity values are negative ({n_rows} of {len(df)} features). "
            "The table looks autoscaled, mean-centred or already log-transformed. PhenoRewire expects "
            "non-negative intensities: it applies its own median normalisation and log2(x + 1), and "
            "treats values <= 0 as not detected. Provide the intensities before scaling."
        )
    elif neg_share > 0:
        chk.warnings.append(
            f"{int(neg.sum())} intensity value(s) are negative; they are treated as not detected."
        )
    elif config.LOG_TRANSFORM and float(np.nanmax(vals)) < LOG_SCALE_MAX:
        chk.warnings.append(
            f"Every intensity is below {LOG_SCALE_MAX:g} (largest: {float(np.nanmax(vals)):.3g}): the table "
            "looks already log-transformed. If it is, set LOG_TRANSFORM: false in the config, or it "
            "will be log-transformed a second time."
        )

    thr = float(config.MIN_NONZERO_INTENSITY)
    Xf = num.fillna(0.0)
    if thr > 0:
        Xf = Xf.where(Xf >= thr, 0.0)
    detected = Xf > 0
    chk.lines.append(f"{_pct(1 - detected.values.mean())} of intensity values are missing or zero")
    presence = detected.mean(axis=1)
    low = presence < float(config.MIN_PRESENCE_GLOBAL)
    if low.all():
        chk.errors.append(
            f"No feature is detected in at least {_pct(config.MIN_PRESENCE_GLOBAL)} of the samples "
            "(MIN_PRESENCE_GLOBAL), so nothing would be left to analyse."
        )
    elif (low & ~all_text).any():
        chk.lines.append(
            f"{int((low & ~all_text).sum())} feature(s) are detected in fewer than "
            f"{_pct(config.MIN_PRESENCE_GLOBAL)} of samples and will be removed (MIN_PRESENCE_GLOBAL)"
        )

    norm = str(config.NORM_METHOD)
    chk.lines.append(
        f"normalisation: {NORM_LABELS[norm]} (NORM_METHOD: {norm})"
        + (", then log2(x + 1)" if config.LOG_TRANSFORM else ", no log transform")
    )
    # Skipped when the table is already refused for negative values: those count as
    # undetected here, and the advice would point at the wrong problem.
    if norm == "median" and not low.all() and neg_share <= NEGATIVE_SHARE_ERROR:
        # Median normalisation counts undetected values: a sample in which most
        # kept features are absent has median 0 and is left unscaled.
        zero_median = int((Xf.loc[~low].median(axis=0) <= 0).sum())
        if zero_median:
            chk.warnings.append(
                f"In {zero_median} sample(s) more than half of the features are not detected, so their "
                "median is 0 and median normalisation cannot scale them. NORM_METHOD: pqn normalises "
                "on the detected features only."
            )
    return cols


def _check_samples(config, meta: pd.DataFrame, cols: list[str], chk: DataCheck) -> Optional[pd.DataFrame]:
    """Matches metadata rows to intensity columns; returns the matched metadata, in column order."""
    sample_col = config.META_SAMPLE_COL
    if sample_col not in meta.columns:
        chk.errors.append(
            f"The metadata has no column '{sample_col}' (META_SAMPLE_COL in the config). "
            f"Its columns are: {_examples(meta.columns, 8)}."
        )
        return None

    meta = meta.copy()
    meta[sample_col] = meta[sample_col].astype(str).str.strip().map(norm_sample_id)
    dup = meta.loc[meta[sample_col].duplicated(keep=False), sample_col].unique().tolist()
    if dup:
        chk.errors.append(f"{len(dup)} sample ID(s) appear more than once in the metadata (e.g. {_examples(dup)}).")
        return None

    in_matrix = matrix_sample_ids(cols, config.INTENSITY_REGEX)
    dup_m = pd.Series(in_matrix)[pd.Series(in_matrix).duplicated(keep=False)].unique().tolist()
    if dup_m:
        chk.errors.append(
            f"{len(dup_m)} sample name(s) occur in more than one intensity column (e.g. {_examples(dup_m)})."
        )
        return None

    meta_ids = set(meta[sample_col])
    missing = [s for s in in_matrix if s not in meta_ids]
    extra = [s for s in meta[sample_col] if s not in set(in_matrix)]
    if len(missing) == len(in_matrix):
        chk.errors.append(
            "No sample of the feature table has a row in the metadata. "
            f"Feature table samples: {_examples(in_matrix, 3)}. "
            f"Metadata '{sample_col}': {_examples(meta[sample_col].tolist(), 3)}."
        )
        return None
    if missing:
        chk.errors.append(
            f"{len(missing)} sample(s) of the feature table have no row in the metadata: {_examples(missing)}."
        )
        return None
    chk.lines.append(f"all {len(in_matrix)} samples of the feature table found in the metadata")
    if extra:
        chk.warnings.append(
            f"{len(extra)} metadata row(s) have no intensity column and are ignored: {_examples(extra)}."
        )
    return meta.set_index(sample_col).loc[in_matrix].reset_index()


def _check_groups(config, aligned: pd.DataFrame, chk: DataCheck) -> Optional[pd.Series]:
    wanted = {k for cond in config.GROUP_DEFINITION.values() for k in cond}
    absent = sorted(k for k in wanted if k not in aligned.columns)
    if absent:
        chk.errors.append(
            f"GROUP_DEFINITION refers to metadata column(s) that do not exist: {_examples(absent)}. "
            f"The metadata columns are: {_examples(aligned.columns, 8)}."
        )
        return None
    groups = pd.Series(assign_groups(aligned, config.GROUP_DEFINITION), index=aligned.index)
    counts = groups.value_counts()
    shown = ", ".join(f"{g} ({int(counts.get(g, 0))})" for g in config.GROUP_DEFINITION)
    chk.lines.append(f"groups: {shown}")
    # assign_groups labels an unmatched sample "Unknown", so test membership rather
    # than nullness: QC injections, blanks and mistyped group values all land here,
    # and a silent drop is exactly what costs a comparison its sample size.
    unmatched = aligned.loc[~groups.isin(list(config.GROUP_DEFINITION)), config.META_SAMPLE_COL]
    if len(unmatched):
        chk.warnings.append(
            f"{len(unmatched)} sample(s) match no group in GROUP_DEFINITION and are not used: "
            f"{_examples(unmatched)}. Expected when the table also holds conditions this "
            f"comparison leaves out, or QC injections and blanks. If it is not, check that the "
            f"values in GROUP_DEFINITION match the metadata exactly."
        )
    if counts.sum() == 0:
        chk.errors.append(
            "No sample matches any group in GROUP_DEFINITION. Check that its values match the metadata "
            "exactly (e.g. 'case' vs 'Case', 1 vs '1')."
        )
        return None
    return groups


def _plan_phenotype(config, groups: pd.Series, chk: DataCheck) -> None:
    mode = config.ANALYSIS_MODE
    with_samples = [g for g in config.GROUP_DEFINITION if int((groups == g).sum()) > 0]
    ref, case = config.PHENO_GROUP_REF, config.PHENO_GROUP_CASE
    explicit = mode in {"phenotype", "both"}
    if ref is not None and case is not None:
        pair = [str(ref), str(case)]
        empty = [g for g in pair if g not in with_samples]
        if empty:
            if explicit or not config.TEMPORAL_CORRELATION.enabled:
                chk.errors.append(
                    f"PHENO_GROUP_REF / PHENO_GROUP_CASE name group(s) with no samples: {_examples(empty)}."
                )
            return
    elif len(with_samples) == 2:
        pair = sorted(with_samples)
    else:
        if explicit:
            chk.errors.append(
                f"A phenotype comparison needs two groups with samples; {len(with_samples)} have samples "
                f"({_examples(with_samples)}). Set PHENO_GROUP_REF and PHENO_GROUP_CASE."
            )
        return
    n_a, n_b = int((groups == pair[0]).sum()), int((groups == pair[1]).sum())
    chk.comparisons.append(
        {"kind": "phenotype", "label": f"{pair[0]} vs {pair[1]}", "states": pair, "n": [n_a, n_b]}
    )


def _plan_temporal(config, aligned: pd.DataFrame, groups: pd.Series, chk: DataCheck) -> None:
    tc = config.TEMPORAL_CORRELATION
    group, time_col = str(tc.phenotype), str(tc.time_column)
    if time_col not in aligned.columns:
        chk.errors.append(
            f"The time column '{time_col}' (TEMPORAL_CORRELATION.time_column) is not in the metadata. "
            f"Its columns are: {_examples(aligned.columns, 8)}."
        )
        return
    in_group = (groups == group).values
    if not in_group.any():
        chk.errors.append(
            f"TEMPORAL_CORRELATION.phenotype = '{group}' has no samples; it must be one of the "
            f"GROUP_DEFINITION names: {_examples(list(config.GROUP_DEFINITION))}."
        )
        return
    times = parse_time_numeric(aligned.loc[in_group, time_col])
    bad = aligned.loc[in_group, time_col][times.isna().values]
    if len(bad):
        chk.errors.append(
            f"{len(bad)} value(s) of '{time_col}' are not numbers (e.g. {_examples(bad.astype(str).unique())}); "
            "write times as numbers, optionally with 'h' (3, 3h)."
        )
        return
    tcounts = times.value_counts().sort_index()
    chk.lines.append(
        f"timepoints in {group}: " + ", ".join(f"{timepoint_label(t)} ({int(n)})" for t, n in tcounts.items())
    )
    t1, t2 = float(config.TIMEPOINT_T1), float(config.TIMEPOINT_T2)
    available = ", ".join(timepoint_label(t) for t in tcounts.index)
    for key, t in (("TIMEPOINT_T1", t1), ("TIMEPOINT_T2", t2)):
        if t not in tcounts.index:
            chk.errors.append(f"{key} = {t:g} has no samples in {group}; available timepoints: {available}.")
    if chk.errors:
        return
    lbl1, lbl2 = timepoint_label(t1), timepoint_label(t2)
    chk.comparisons.append(
        {
            "kind": "temporal",
            "label": f"{group}: {lbl1} vs {lbl2}",
            "states": [lbl1, lbl2],
            "n": [int(tcounts[t1]), int(tcounts[t2])],
            "group": group,
        }
    )


def check_inputs(
    config,
    df: Optional[pd.DataFrame] = None,
    meta: Optional[pd.DataFrame] = None,
) -> DataCheck:
    """Check the input tables named in ``config``; never raises on bad data."""
    chk = DataCheck(data_matrix=_shown(config.DATA_MATRIX), metadata=_shown(config.METADATA))
    if df is None:
        df = pd.read_csv(config.DATA_MATRIX, low_memory=False)
        df.columns = [str(c).strip() for c in df.columns]
    if meta is None:
        meta = pd.read_csv(config.METADATA, low_memory=False)
        meta.columns = [str(c).strip() for c in meta.columns]

    cols = _check_feature_table(config, df, chk)
    if cols is None:
        return chk
    aligned = _check_samples(config, meta, cols, chk)
    if aligned is None:
        return chk
    groups = _check_groups(config, aligned, chk)
    if groups is None:
        return chk

    mode = config.ANALYSIS_MODE
    if mode in {"auto", "phenotype", "both"}:
        _plan_phenotype(config, groups, chk)
    if config.TEMPORAL_CORRELATION.enabled and mode in {"auto", "temporal", "both"}:
        _plan_temporal(config, aligned, groups, chk)
    if not chk.comparisons and not chk.errors:
        chk.errors.append(
            "Nothing to compare: no two groups with samples for a phenotype comparison, and the "
            "temporal analysis is not enabled."
        )

    for comp in chk.comparisons:
        op = operating_range(*comp["n"])
        comp["operating_range"] = op
        sizes = f"{comp['n'][0]} vs {comp['n'][1]} samples"
        if op["status"] == "within":
            chk.lines.append(f"comparison {comp['label']}: {sizes}, within the operating range")
        elif op["status"] == "reduced":
            chk.lines.append(f"comparison {comp['label']}: {sizes}, reduced range")
            chk.warnings.append(
                f"{comp['label']}: {op['min_n_per_state']} samples in the smaller state; "
                f"{OPERATING_RANGE_RECOMMENDED_N} or more are recommended. Expect low recall."
            )
        else:
            chk.lines.append(f"comparison {comp['label']}: {sizes}, below the operating range")
            chk.warnings.append(
                f"{comp['label']}: {op['min_n_per_state']} samples in the smaller state, below the "
                f"operating range ({OPERATING_RANGE_MIN_N}). The run will continue, but rewiring at this "
                "size is not reliably distinguishable from chance."
            )
    return chk


def _shown(path) -> str:
    path = Path(path)
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve())).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")
