#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
CUR = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
PROT = BASE / "results" / "features" / "protparam_features.csv"
STR = BASE / "results" / "features" / "structural_features_from_models.csv"
MOD = BASE / "models" / "swiss_model" / "swiss_model_manifest_registered.csv"
OUT = BASE / "results" / "ml" / "ml_feature_matrix.csv"

cur = pd.read_csv(CUR)
prot = pd.read_csv(PROT)
st = pd.read_csv(STR)
mod = pd.read_csv(MOD)

st["uniprot_accession"] = st["pdb"].astype(str).str.replace("_swissmodel.pdb", "", regex=False)
mod = mod[["uniprot_accession","swiss_model_status","model_path"]]

merged = cur.merge(prot, on=["uniprot_accession","organism","organism_type","gh_family","length","optimum_temperature","optimum_pH"], how="left")
merged = merged.merge(st[["uniprot_accession","hbond_proxy","salt_bridge_proxy","disulfide_count"]], on="uniprot_accession", how="left")
merged = merged.merge(mod, on="uniprot_accession", how="left")

# scaffold labels for current project
merged["thermostability_label"] = merged["optimum_temperature"].fillna("").astype(str).str.strip()
merged["cellulolytic_activity_label"] = ""

OUT.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(merged)}")
