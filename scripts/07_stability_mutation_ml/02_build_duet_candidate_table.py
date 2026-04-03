#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
CUR = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
STRUCT = BASE / "results" / "features" / "structural_features_from_models.csv"
PROT = BASE / "results" / "features" / "protparam_features.csv"
MODELS = BASE / "models" / "swiss_model" / "swiss_model_manifest_registered.csv"
OUT = BASE / "results" / "mutations" / "duet_candidate_table.csv"

# load
cur = pd.read_csv(CUR).copy()
sdf = pd.read_csv(STRUCT).copy()
pdf = pd.read_csv(PROT).copy()
mdf = pd.read_csv(MODELS).copy()

# normalize ids
cur["uniprot_accession"] = cur["uniprot_accession"].fillna("").astype(str).str.strip()
pdf["uniprot_accession"] = pdf["uniprot_accession"].fillna("").astype(str).str.strip()
mdf["uniprot_accession"] = mdf["uniprot_accession"].fillna("").astype(str).str.strip()

sdf["uniprot_accession"] = (
    sdf["pdb"].astype(str)
    .str.replace("_swissmodel.pdb", "", regex=False)
    .str.strip()
)

# keep only received models
mdf = mdf[mdf["swiss_model_status"] == "model_received"].copy()

# use curated table as the base table so organism/org_type/gh remain present
base = cur[[
    "uniprot_accession",
    "organism",
    "organism_type",
    "gh_family",
    "optimum_temperature",
    "optimum_pH",
]].copy()

# reduce other tables to needed columns
pdf_keep = [c for c in [
    "uniprot_accession",
    "gravy",
    "predicted_pI",
    "instability_index",
    "aromaticity",
    "molecular_weight",
] if c in pdf.columns]

sdf_keep = [c for c in [
    "uniprot_accession",
    "hbond_proxy",
    "salt_bridge_proxy",
    "disulfide_count",
] if c in sdf.columns]

mdf_keep = [c for c in [
    "uniprot_accession",
    "model_path",
    "template",
    "gmqe",
    "qmean",
] if c in mdf.columns]

merged = base.merge(pdf[pdf_keep], on="uniprot_accession", how="left")
merged = merged.merge(sdf[sdf_keep], on="uniprot_accession", how="left")
merged = merged.merge(mdf[mdf_keep], on="uniprot_accession", how="inner")

# final DUET-ready table
out = merged.copy()
out["mutation_chain"] = "A"
out["mutation"] = ""
out["mutation_rationale"] = ""
out["duet_ddg"] = ""

OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(out)}")
print(out[[
    "uniprot_accession",
    "organism_type",
    "gh_family",
    "model_path"
]].head(10).to_string(index=False))
