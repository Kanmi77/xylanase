#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
MANIFEST = BASE / "models" / "swiss_model" / "swiss_model_manifest_final.csv"
OUTDIR = BASE / "models" / "swiss_model" / "outputs"
OUTCSV = BASE / "models" / "swiss_model" / "swiss_model_manifest_registered.csv"

manifest = pd.read_csv(MANIFEST)
OUTDIR.mkdir(parents=True, exist_ok=True)

# Always rebuild model_path from actual files present in outputs
manifest["model_path"] = ""

pdb_map = {}
for pdb_file in OUTDIR.glob("*.pdb"):
    acc = pdb_file.stem.split("_")[0].strip()
    pdb_map[acc] = str(pdb_file)

manifest["model_path"] = manifest["uniprot_accession"].astype(str).str.strip().map(pdb_map).fillna("")
manifest["swiss_model_status"] = manifest["model_path"].apply(
    lambda x: "model_received" if str(x).strip() else "pending_submission"
)

OUTCSV.parent.mkdir(parents=True, exist_ok=True)
manifest.to_csv(OUTCSV, index=False)

print(f"Saved: {OUTCSV}")
print(manifest["swiss_model_status"].value_counts(dropna=False))
print(f"PDB files found in outputs: {len(pdb_map)}")
