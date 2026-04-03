#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
MODELS = BASE / "models" / "swiss_model" / "swiss_model_manifest_registered.csv"
OUT = BASE / "results" / "docking" / "docking_manifest.csv"

df = pd.read_csv(MODELS)
df = df[df["swiss_model_status"] == "model_received"].copy()

out = df[["uniprot_accession","organism","organism_type","gh_family","model_path"]].copy()
out["ligand"] = "cellobiose_or_cellotriose"
out["vina_score"] = ""
out["notes"] = ""

OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(out)}")
