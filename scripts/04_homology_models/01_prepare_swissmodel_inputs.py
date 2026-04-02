#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

base = Path.home() / "xylanase"
master_file = base / "data" / "curated" / "xylanase_master_table.csv"
outdir = base / "models" / "swiss_model" / "inputs"
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(master_file)

no_pdb = df[df["has_pdb"] == 0].copy()

fasta_file = outdir / "sequences_without_pdb.fasta"
with open(fasta_file, "w", encoding="utf-8") as f:
    for _, row in no_pdb.iterrows():
        acc = str(row["uniprot_accession"])
        org = str(row["organism"]).replace(" ", "_")
        seq = str(row["sequence"])
        f.write(f">{acc}|{org}\n{seq}\n")

manifest = outdir / "swiss_model_manifest.csv"
no_pdb.to_csv(manifest, index=False)

print(f"Saved FASTA for SWISS-MODEL: {fasta_file}")
print(f"Saved manifest: {manifest}")
print(f"Sequences lacking PDB: {len(no_pdb)}")
