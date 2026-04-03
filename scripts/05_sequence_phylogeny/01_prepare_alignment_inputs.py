#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
OUTDIR = BASE / "data" / "phylogeny"

OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE).copy()

df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()
df["organism_type"] = df["organism_type"].fillna("").astype(str).str.strip().str.lower()
df["gh_family"] = df["gh_family"].fillna("").astype(str).str.strip()
df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()
df["optimum_temperature"] = df["optimum_temperature"].fillna("").astype(str).str.strip()

# simple temperature class for labeled entries only
def classify_temp(x):
    if not x:
        return "unlabeled"
    try:
        val = float(x)
    except Exception:
        return "unlabeled"
    if val >= 60:
        return "thermostable"
    return "mesophilic"

df["temp_class"] = df["optimum_temperature"].apply(classify_temp)

# save metadata table for downstream plotting
meta_cols = [
    "uniprot_accession",
    "organism",
    "organism_type",
    "gh_family",
    "optimum_temperature",
    "optimum_pH",
    "temp_class",
]
df[meta_cols].to_csv(OUTDIR / "phylogeny_metadata.csv", index=False)

# global fasta
with open(OUTDIR / "all_filtered_xylanases.fasta", "w", encoding="utf-8") as f:
    for _, row in df.iterrows():
        if row["uniprot_accession"] and row["sequence"]:
            header = f'>{row["uniprot_accession"]}|{row["organism_type"]}|{row["gh_family"]}|{row["temp_class"]}'
            f.write(header + "\n")
            f.write(row["sequence"] + "\n")

# subgroup fastas
groups = [
    ("bacterial", "GH10"),
    ("bacterial", "GH11"),
    ("fungal", "GH10"),
    ("fungal", "GH11"),
]

for org, gh in groups:
    sub = df[(df["organism_type"] == org) & (df["gh_family"] == gh)].copy()
    out = OUTDIR / f"{org}_{gh}.fasta"
    with open(out, "w", encoding="utf-8") as f:
        for _, row in sub.iterrows():
            if row["uniprot_accession"] and row["sequence"]:
                header = f'>{row["uniprot_accession"]}|{row["organism_type"]}|{row["gh_family"]}|{row["temp_class"]}'
                f.write(header + "\n")
                f.write(row["sequence"] + "\n")

print("Saved FASTA and metadata files to:", OUTDIR)
print("Total rows:", len(df))
print(df[["organism_type", "gh_family"]].value_counts())
print("\nTemperature class counts:")
print(df["temp_class"].value_counts(dropna=False))
