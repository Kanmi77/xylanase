#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

base = Path("/home/ubuntu/xylanase/xylanase")

infile = base / "data/curated/xylanase_master_table.csv"
outfile = base / "data/curated/xylanase_master_table_annotated.csv"
report = base / "results/reports/curation_results_after_cazy.txt"

df = pd.read_csv(infile)

def infer_organism_type(tax):
    tax = str(tax).lower()
    if "bacteria" in tax:
        return "bacterial"
    if "fungi" in tax:
        return "fungal"
    return np.nan

def infer_gh(row):
    text = " ".join([
        str(row.get("protein_name", "")),
        str(row.get("entry_name", "")),
        str(row.get("gene_names", "")),
    ]).upper()

    if "GH10" in text or "XYN10" in text or "ENDO-1,4-BETA-XYLANASE GH10" in text:
        return "GH10"
    if "GH11" in text or "XYN11" in text or "ENDO-1,4-BETA-XYLANASE GH11" in text:
        return "GH11"
    return np.nan

if "organism_type" not in df.columns:
    df["organism_type"] = np.nan
if "gh_family" not in df.columns:
    df["gh_family"] = np.nan

df["organism_type"] = df["organism_type"].fillna(df["taxonomy"].apply(infer_organism_type))
df["gh_family"] = df.apply(
    lambda row: row["gh_family"] if pd.notna(row["gh_family"]) else infer_gh(row),
    axis=1
)

df.to_csv(outfile, index=False)

with open(report, "w") as fh:
    fh.write("Curation results after GH/name/taxonomy annotation\n\n")
    fh.write("GH family counts:\n")
    fh.write(df["gh_family"].value_counts(dropna=False).to_string())
    fh.write("\n\nOrganism type counts:\n")
    fh.write(df["organism_type"].value_counts(dropna=False).to_string())
    fh.write("\n")

print(f"Saved: {outfile}")
print(f"Saved: {report}")
print("\nGH family:")
print(df["gh_family"].value_counts(dropna=False))
print("\nOrganism type:")
print(df["organism_type"].value_counts(dropna=False))
