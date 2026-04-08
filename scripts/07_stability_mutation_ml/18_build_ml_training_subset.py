#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "ml" / "ml_feature_matrix_with_thermoprot.csv"
OUTFILE = BASE / "results" / "ml" / "ml_training_subset.csv"

df = pd.read_csv(INFILE).copy()

feature_cols = [c for c in [
    "length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "gravy",
    "predicted_pI",
    "secondary_structure_fraction_helix",
    "secondary_structure_fraction_turn",
    "secondary_structure_fraction_sheet",
    "hbond_proxy",
    "salt_bridge_proxy",
    "disulfide_count",
    "Probability"
] if c in df.columns]

keep_cols = ["uniprot_accession", "organism_type", "gh_family"] + feature_cols

out = df[keep_cols].copy()
out["cellulolytic_activity_label"] = ""

out.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(out)}")
print("Feature columns used:")
print(feature_cols)
