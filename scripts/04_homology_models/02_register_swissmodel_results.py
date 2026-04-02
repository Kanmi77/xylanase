#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
MANIFEST = BASE / "models" / "swiss_model" / "swiss_model_manifest_final.csv"
OUTDIR = BASE / "models" / "swiss_model" / "outputs"
OUTCSV = BASE / "models" / "swiss_model" / "swiss_model_manifest_registered.csv"

manifest = pd.read_csv(MANIFEST)
OUTDIR.mkdir(parents=True, exist_ok=True)

pdb_map = {}
for pdb_file in OUTDIR.glob("*.pdb"):
    acc = pdb_file.stem.split("_")[0]
    pdb_map[acc] = str(pdb_file)

manifest["model_path"] = manifest["uniprot_accession"].map(pdb_map).fillna(manifest["model_path"])
manifest["swiss_model_status"] = manifest["model_path"].apply(
    lambda x: "model_received" if str(x).strip() else "pending_submission"
)

manifest.to_csv(OUTCSV, index=False)
print(f"Saved: {OUTCSV}")
print(manifest["swiss_model_status"].value_counts(dropna=False))
