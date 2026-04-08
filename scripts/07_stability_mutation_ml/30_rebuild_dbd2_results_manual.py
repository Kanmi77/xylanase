#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "dbd2_input_list.csv"
OUTFILE = BASE / "results" / "mutations" / "dbd2_results_manual.csv"

df = pd.read_csv(INFILE).copy()

# ensure manual result columns exist
df["dbd2_pair"] = df.get("dbd2_pair", "")
df["dbd2_score"] = ""
df["dbd2_prediction"] = ""
df["dbd2_interpretation"] = ""

# keep clean order
cols = [c for c in [
    "uniprot_accession",
    "organism",
    "organism_type",
    "gh_family",
    "model_path",
    "dbd2_pair",
    "dbd2_score",
    "dbd2_prediction",
    "dbd2_interpretation"
] if c in df.columns]

df = df[cols].copy()
df.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(df)}")
print(df.head(10).to_string(index=False))
