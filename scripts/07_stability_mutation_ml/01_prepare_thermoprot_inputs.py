#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
OUTFA = BASE / "results" / "thermostability" / "thermoprot_input.fasta"
META = BASE / "results" / "thermostability" / "thermoprot_input_metadata.csv"

df = pd.read_csv(INFILE).copy()
df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()
df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()

OUTFA.parent.mkdir(parents=True, exist_ok=True)
with open(OUTFA, "w", encoding="utf-8") as f:
    for _, row in df.iterrows():
        acc = row["uniprot_accession"]
        seq = row["sequence"]
        if acc and seq:
            f.write(f">{acc}\n{seq}\n")

df[["uniprot_accession","organism","organism_type","gh_family","optimum_temperature","optimum_pH"]].to_csv(META, index=False)

print(f"Saved: {OUTFA}")
print(f"Saved: {META}")
print(f"Rows: {len(df)}")
