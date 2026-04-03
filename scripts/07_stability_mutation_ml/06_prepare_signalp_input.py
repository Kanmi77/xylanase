#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
OUTFA = BASE / "results" / "ml" / "signalp_input.fasta"

df = pd.read_csv(INFILE).copy()
df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()
df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()

with open(OUTFA, "w", encoding="utf-8") as f:
    for _, row in df.iterrows():
        acc = row["uniprot_accession"]
        seq = row["sequence"]
        if acc and seq:
            f.write(f">{acc}\n{seq}\n")

print(f"Saved: {OUTFA}")
