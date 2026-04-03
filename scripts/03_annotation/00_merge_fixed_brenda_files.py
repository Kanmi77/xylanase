#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase" / "data" / "brenda" / "processed"

files = [
    BASE / "xylanase_temperature_optimum_FIXED.csv",
    BASE / "xylanase_temperature_range_FIXED.csv",
    BASE / "xylanase_temperature_stability_FIXED.csv",
]

dfs = []
for f in files:
    df = pd.read_csv(f)
    if "uniprot" in df.columns and "uniprot_accession" not in df.columns:
        df = df.rename(columns={"uniprot": "uniprot_accession"})
    dfs.append(df)

merged = pd.concat(dfs, ignore_index=True)
out = BASE / "brenda_xylanase_structured.csv"
merged.to_csv(out, index=False)

print("Saved:", out)
print("Rows:", len(merged))
print("Columns:", merged.columns.tolist())
