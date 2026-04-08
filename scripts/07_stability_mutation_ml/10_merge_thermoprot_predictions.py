#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
ML = BASE / "results" / "ml" / "ml_feature_matrix.csv"
PRED = BASE / "results" / "thermostability" / "predictions" / "thermoprot_predictions.csv"
OUT = BASE / "results" / "ml" / "ml_feature_matrix_with_thermoprot.csv"

ml = pd.read_csv(ML)
pred = pd.read_csv(PRED)

ml["uniprot_accession"] = ml["uniprot_accession"].astype(str).str.strip()

# ThermoProt uses FASTA header output
pred["Header"] = pred["Header"].astype(str).str.strip()

pred["uniprot_accession"] = (
    pred["Header"]
    .str.replace("^>", "", regex=True)
    .str.split().str[0]
    .str.split("|").str[0]
    .str.strip()
)

out = ml.merge(pred, on="uniprot_accession", how="left")
out.to_csv(OUT, index=False)

print("ThermoProt ID column used: Header")
print(f"Saved: {OUT}")
print(f"Rows: {len(out)}")
