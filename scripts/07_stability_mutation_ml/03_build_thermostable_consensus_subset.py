#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
OUTFA = BASE / "results" / "mutations" / "thermostable_labeled_subset.fasta"
OUTCSV = BASE / "results" / "mutations" / "thermostable_labeled_subset.csv"

df = pd.read_csv(INFILE).copy()
df["optimum_temperature"] = df["optimum_temperature"].fillna("").astype(str).str.strip()
df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()
df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()

def temp_class(x):
    try:
        return "thermostable" if float(x) >= 60 else "mesophilic"
    except Exception:
        return "unlabeled"

df["temp_class"] = df["optimum_temperature"].apply(temp_class)
sub = df[df["temp_class"] == "thermostable"].copy()

sub.to_csv(OUTCSV, index=False)

with open(OUTFA, "w", encoding="utf-8") as f:
    for _, row in sub.iterrows():
        if row["uniprot_accession"] and row["sequence"]:
            f.write(f'>{row["uniprot_accession"]}\n{row["sequence"]}\n')

print(f"Saved: {OUTCSV}")
print(f"Saved: {OUTFA}")
print(f"Thermostable labeled rows: {len(sub)}")
