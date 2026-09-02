"""
Generate synthetic example data for PhenoRewire.

Run once to create feature_table.csv and metadata.csv in this directory:
    python example/generate_example_data.py

The script produces a small but realistic metabolomics dataset with:
  - 22 features (20 biologically patterned + 2 pure noise)
  - 20 samples (10 reference, 10 case)
  - Clear correlation structure that produces visible rewiring
"""

import numpy as np
import pandas as pd
from pathlib import Path

rng = np.random.default_rng(42)
OUT = Path(__file__).parent

N_REF = 10
N_CASE = 10
N = N_REF + N_CASE

sample_ids = [f"Sample_{i:02d}" for i in range(1, N + 1)]
group_labels = ["reference"] * N_REF + ["case"] * N_CASE

# ── Reference latent factors ──────────────────────────────────────────────────
# Two independent metabolic modules active in reference
ref_mod1 = rng.standard_normal(N_REF)   # Module A (energy metabolism)
ref_mod2 = rng.standard_normal(N_REF)   # Module B (amino acid metabolism)

# ── Case latent factors ───────────────────────────────────────────────────────
# Module A is attenuated; Module C replaces Module B; new Module D emerges
case_mod1 = rng.standard_normal(N_CASE) * 0.4   # Module A dampened in case
case_mod3 = rng.standard_normal(N_CASE)          # Module C (case-specific)
case_mod4 = rng.standard_normal(N_CASE)          # Module D (case-specific)

def noise(n):
    return rng.normal(0, 0.25, n)

# ── Feature intensities ───────────────────────────────────────────────────────
# F001–F004: Module A — correlated in both groups, slightly higher in case
f001_ref = 8.5 + 1.2 * ref_mod1 + noise(N_REF);  f001_case = 9.8 + 0.5 * case_mod1 + noise(N_CASE)
f002_ref = 7.2 + 1.0 * ref_mod1 + noise(N_REF);  f002_case = 8.4 + 0.4 * case_mod1 + noise(N_CASE)
f003_ref = 6.8 + 0.9 * ref_mod1 + noise(N_REF);  f003_case = 7.9 + 0.4 * case_mod1 + noise(N_CASE)
f004_ref = 9.1 + 1.1 * ref_mod1 + noise(N_REF);  f004_case = 10.2 + 0.5 * case_mod1 + noise(N_CASE)

# F005–F008: Module B in reference, not correlated in case (rewired: ref-specific edges)
f005_ref = 7.5 + 1.3 * ref_mod2 + noise(N_REF);  f005_case = 8.1 + noise(N_CASE)
f006_ref = 6.9 + 1.1 * ref_mod2 + noise(N_REF);  f006_case = 7.4 + noise(N_CASE)
f007_ref = 8.0 + 1.0 * ref_mod2 + noise(N_REF);  f007_case = 8.6 + noise(N_CASE)
f008_ref = 7.3 + 0.9 * ref_mod2 + noise(N_REF);  f008_case = 8.0 + noise(N_CASE)

# F009–F012: Module C in case, not correlated in reference (rewired: case-specific edges)
f009_ref = 5.5 + noise(N_REF);  f009_case = 6.8 + 1.2 * case_mod3 + noise(N_CASE)
f010_ref = 6.1 + noise(N_REF);  f010_case = 7.5 + 1.0 * case_mod3 + noise(N_CASE)
f011_ref = 5.8 + noise(N_REF);  f011_case = 7.1 + 1.1 * case_mod3 + noise(N_CASE)
f012_ref = 6.3 + noise(N_REF);  f012_case = 7.8 + 0.9 * case_mod3 + noise(N_CASE)

# F013–F015: Module D in case, higher than reference (strongly differentially abundant)
f013_ref = 4.2 + noise(N_REF);  f013_case = 7.5 + 1.3 * case_mod4 + noise(N_CASE)
f014_ref = 3.9 + noise(N_REF);  f014_case = 7.1 + 1.1 * case_mod4 + noise(N_CASE)
f015_ref = 4.5 + noise(N_REF);  f015_case = 7.8 + 1.2 * case_mod4 + noise(N_CASE)

# F016–F017: Sign-switch — positively correlated in reference, negatively in case
f016_ref = 7.0 + 1.0 * ref_mod1 + noise(N_REF);  f016_case = 8.5 - 0.9 * case_mod1 + noise(N_CASE)
f017_ref = 6.5 + 0.9 * ref_mod1 + noise(N_REF);  f017_case = 7.9 - 0.8 * case_mod1 + noise(N_CASE)

# F018–F019: Moderately associated with phenotype (borderline)
f018_ref = 5.0 + noise(N_REF);  f018_case = 6.2 + noise(N_CASE)
f019_ref = 4.8 + noise(N_REF);  f019_case = 5.7 + noise(N_CASE)

# F020–F022: Pure noise (not phenotype-associated, should not be selected)
f020 = 6.0 + rng.standard_normal(N)
f021 = 5.5 + rng.standard_normal(N)
f022 = 7.0 + rng.standard_normal(N)

# ── Assemble feature table ────────────────────────────────────────────────────
intensities = {
    "F001": np.concatenate([f001_ref, f001_case]),
    "F002": np.concatenate([f002_ref, f002_case]),
    "F003": np.concatenate([f003_ref, f003_case]),
    "F004": np.concatenate([f004_ref, f004_case]),
    "F005": np.concatenate([f005_ref, f005_case]),
    "F006": np.concatenate([f006_ref, f006_case]),
    "F007": np.concatenate([f007_ref, f007_case]),
    "F008": np.concatenate([f008_ref, f008_case]),
    "F009": np.concatenate([f009_ref, f009_case]),
    "F010": np.concatenate([f010_ref, f010_case]),
    "F011": np.concatenate([f011_ref, f011_case]),
    "F012": np.concatenate([f012_ref, f012_case]),
    "F013": np.concatenate([f013_ref, f013_case]),
    "F014": np.concatenate([f014_ref, f014_case]),
    "F015": np.concatenate([f015_ref, f015_case]),
    "F016": np.concatenate([f016_ref, f016_case]),
    "F017": np.concatenate([f017_ref, f017_case]),
    "F018": np.concatenate([f018_ref, f018_case]),
    "F019": np.concatenate([f019_ref, f019_case]),
    "F020": f020,
    "F021": f021,
    "F022": f022,
}

annotations = [
    ("F001", 133.014, 0.98,  "lactate"),
    ("F002", 175.119, 2.31,  "citrulline"),
    ("F003", 146.058, 1.87,  "glutamate"),
    ("F004", 132.054, 3.12,  "asparagine"),
    ("F005", 204.090, 4.55,  "adenosine"),
    ("F006", 267.097, 5.21,  "AMP"),
    ("F007", 348.069, 4.83,  "ADP"),
    ("F008", 428.037, 4.91,  "ATP"),
    ("F009", 117.043, 1.44,  "threonine"),
    ("F010", 131.058, 1.69,  "methionine"),
    ("F011", 147.053, 1.55,  "glutamine"),
    ("F012", 132.102, 2.04,  "leucine"),
    ("F013", 168.065, 3.77,  "phenylalanine"),
    ("F014", 182.081, 4.22,  "tyrosine"),
    ("F015", 204.090, 5.01,  "tryptophan"),
    ("F016", 191.020, 1.11,  "citrate"),
    ("F017", 173.009, 0.95,  "isocitrate"),
    ("F018", 90.032,  0.72,  "pyruvate"),
    ("F019", 116.047, 1.33,  "succinate"),
    ("F020", 245.077, 6.14,  None),
    ("F021", 309.045, 7.88,  None),
    ("F022", 181.061, 8.33,  None),
]

rows = []
for feat_id, mz, rt, name in annotations:
    row = {"feature_id": feat_id, "mz": mz, "rt": rt, "name": name if name else ""}
    for sid, val in zip(sample_ids, intensities[feat_id]):
        row[sid] = round(max(0.0, float(val)), 4)
    rows.append(row)

feature_table = pd.DataFrame(rows)
feature_table.to_csv(OUT / "feature_table.csv", index=False)
print(f"Written: {OUT / 'feature_table.csv'}  ({feature_table.shape[0]} features × {N} samples)")

# ── Metadata ──────────────────────────────────────────────────────────────────
metadata = pd.DataFrame({
    "sample_id": sample_ids,
    "group_label": group_labels,
})
metadata.to_csv(OUT / "metadata.csv", index=False)
print(f"Written: {OUT / 'metadata.csv'}  ({metadata.shape[0]} samples)")
