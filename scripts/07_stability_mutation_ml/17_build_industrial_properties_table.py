#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
CUR = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
PROT = BASE / "results" / "features" / "protparam_features.csv"
THERMO = BASE / "results" / "thermostability" / "predictions" / "thermoprot_predictions.csv"
OUT = BASE / "results" / "reports" / "industrial_properties_table.csv"

cur = pd.read_csv(CUR)
prot = pd.read_csv(PROT)
thermo = pd.read_csv(THERMO)

cur["uniprot_accession"] = cur["uniprot_accession"].astype(str).str.strip()
prot["uniprot_accession"] = prot["uniprot_accession"].astype(str).str.strip()

thermo["Header"] = thermo["Header"].astype(str).str.strip()
thermo["uniprot_accession"] = (
    thermo["Header"]
    .str.replace("^>", "", regex=True)
    .str.split().str[0]
    .str.split("|").str[0]
    .str.strip()
)

merged = cur.merge(
    prot[["uniprot_accession", "predicted_pI", "gravy", "instability_index"]],
    on="uniprot_accession",
    how="left"
).merge(
    thermo[["uniprot_accession", "Probability", "Class", "Prediction"]],
    on="uniprot_accession",
    how="left"
)

out = merged[[
    "uniprot_accession",
    "organism",
    "organism_type",
    "gh_family",
    "optimum_temperature",
    "optimum_pH",
    "predicted_pI",
    "gravy",
    "instability_index",
    "Probability",
    "Class",
    "Prediction"
]].copy()

out["ionic_conditions"] = ""
out["secretion_signal"] = ""

OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(out)}")
