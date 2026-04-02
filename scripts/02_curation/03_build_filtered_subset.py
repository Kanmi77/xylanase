#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_master_table_annotated.csv"
OUTFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11.csv"

df = pd.read_csv(INFILE)

df["gh_family"] = df["gh_family"].fillna("").astype(str).str.strip()
df["organism_type"] = df["organism_type"].fillna("").astype(str).str.strip().str.lower()

filtered = df[
    df["gh_family"].isin(["GH10", "GH11"]) &
    df["organism_type"].isin(["bacterial", "fungal"])
].copy()

filtered.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows in filtered subset: {len(filtered)}")
print("\nBy GH family:")
print(filtered["gh_family"].value_counts(dropna=False))
print("\nBy organism type:")
print(filtered["organism_type"].value_counts(dropna=False))
print("\nBy organism_type x gh_family:")
print(pd.crosstab(filtered["organism_type"], filtered["gh_family"]))
