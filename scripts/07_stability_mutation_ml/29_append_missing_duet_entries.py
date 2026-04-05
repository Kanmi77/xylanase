#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

OUT = Path("/home/ubuntu/xylanase/xylanase/results/mutations/duet_results_manual.csv")

rows = [
    {
        "uniprot_accession": "A0A5P6VQ80",
        "header_mutation": "G317V",
        "mutation": "G317V",
        "tested_mutation": "G317V",
        "mutation_match": "match",
        "original_aa": "GLY",
        "position": 317,
        "new_aa": "VAL",
        "chain": "A",
        "secondary_structure": "Loop or irregular",
        "mcsm_ddg": -0.637,
        "mcsm_interpretation": "Destabilizing",
        "sdm_ddg": 0.52,
        "sdm_interpretation": "Stabilizing",
        "duet_ddg": -0.211,
        "duet_interpretation": "Destabilizing",
        "raw_block": """A0A5P6VQ80 - G317V
mCSM Predicted Stability Change (ΔΔG):
-0.637 kcal/mol (Destabilizing)

SDM Predicted Stability Change (ΔΔG):
0.52 kcal/mol (Stabilizing)

DUET Predicted Stability Change (ΔΔG):
-0.211 kcal/mol (Destabilizing)

Mutation:
Wild-type: GLY
Position: 317
Mutant-type: VAL
Chain: A
Secondary structure: Loop or irregular"""
    },
    {
        "uniprot_accession": "A0A286N5W6",
        "header_mutation": "S445D",
        "mutation": "S445D",
        "tested_mutation": "S445D",
        "mutation_match": "match",
        "original_aa": "SER",
        "position": 445,
        "new_aa": "ASP",
        "chain": "A",
        "secondary_structure": "Loop or irregular",
        "mcsm_ddg": -1.275,
        "mcsm_interpretation": "Destabilizing",
        "sdm_ddg": -0.05,
        "sdm_interpretation": "Destabilizing",
        "duet_ddg": -1.011,
        "duet_interpretation": "Destabilizing",
        "raw_block": """A0A286N5W6 - S445D
mCSM Predicted Stability Change (ΔΔG):
-1.275 kcal/mol (Destabilizing)

SDM Predicted Stability Change (ΔΔG):
-0.05 kcal/mol (Destabilizing)

DUET Predicted Stability Change (ΔΔG):
-1.011 kcal/mol (Destabilizing)

Mutation:
Wild-type: SER
Position: 445
Mutant-type: ASP
Chain: A
Secondary structure: Loop or irregular"""
    },
    {
        "uniprot_accession": "E0S2F5",
        "header_mutation": "D453S",
        "mutation": "D453S",
        "tested_mutation": "D453S",
        "mutation_match": "match",
        "original_aa": "ASP",
        "position": 453,
        "new_aa": "SER",
        "chain": "A",
        "secondary_structure": "Loop or irregular",
        "mcsm_ddg": -1.405,
        "mcsm_interpretation": "Destabilizing",
        "sdm_ddg": -1.08,
        "sdm_interpretation": "Destabilizing",
        "duet_ddg": -1.421,
        "duet_interpretation": "Destabilizing",
        "raw_block": """E0S2F5 - D453S
mCSM Predicted Stability Change (ΔΔG):
-1.405 kcal/mol (Destabilizing)

SDM Predicted Stability Change (ΔΔG):
-1.08 kcal/mol (Destabilizing)

DUET Predicted Stability Change (ΔΔG):
-1.421 kcal/mol (Destabilizing)

Mutation:
Wild-type: ASP
Position: 453
Mutant-type: SER
Chain: A
Secondary structure: Loop or irregular"""
    },
]

df = pd.read_csv(OUT)
add = pd.DataFrame(rows)

# align columns
for col in add.columns:
    if col not in df.columns:
        df[col] = pd.NA
for col in df.columns:
    if col not in add.columns:
        add[col] = pd.NA

combined = pd.concat([df[df.columns], add[df.columns]], ignore_index=True)
combined = combined.drop_duplicates(subset=["uniprot_accession", "mutation"], keep="first")
combined = combined.sort_values(["uniprot_accession", "position"], na_position="last").reset_index(drop=True)

combined.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print("Rows:", len(combined))
print("\nPer protein:")
print(combined["uniprot_accession"].value_counts().to_string())
