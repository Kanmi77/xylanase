#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"

DUET = BASE / "results" / "mutations" / "duet_results_manual.csv"
DBD2 = BASE / "results" / "mutations" / "dbd2_results_manual.csv"
FOCUSED = BASE / "results" / "mutations" / "focused_mutation_candidates.csv"
OUT = BASE / "results" / "mutations" / "final_mutation_results_integrated.csv"

focused = pd.read_csv(FOCUSED).copy()
duet = pd.read_csv(DUET).copy() if DUET.exists() else pd.DataFrame()
dbd2 = pd.read_csv(DBD2).copy() if DBD2.exists() else pd.DataFrame()

for df in [focused, duet, dbd2]:
    if not df.empty and "uniprot_accession" in df.columns:
        df["uniprot_accession"] = df["uniprot_accession"].astype(str).str.strip()

# merge DUET if present
if not duet.empty:
    duet_keep = [c for c in [
        "uniprot_accession",
        "tested_mutation",
        "duet_ddg",
        "duet_interpretation"
    ] if c in duet.columns]
    duet = duet[duet_keep].copy()
    focused = focused.merge(duet, on="uniprot_accession", how="left")

# merge DbD2 if present
if not dbd2.empty:
    dbd2_keep = [c for c in [
        "uniprot_accession",
        "dbd2_pair",
        "dbd2_score",
        "dbd2_prediction",
        "dbd2_interpretation"
    ] if c in dbd2.columns]
    dbd2 = dbd2[dbd2_keep].copy()
    focused = focused.merge(dbd2, on="uniprot_accession", how="left")

OUT.parent.mkdir(parents=True, exist_ok=True)
focused.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(focused)}")
print(focused.head(20).to_string(index=False))
