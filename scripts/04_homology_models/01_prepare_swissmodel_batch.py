#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_annotated.csv"
if not INFILE.exists():
    INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11.csv"

OUTDIR = BASE / "models" / "swiss_model" / "inputs"
MANIFEST = BASE / "models" / "swiss_model" / "swiss_model_manifest_final.csv"

df = pd.read_csv(INFILE)

df["pdb_ids"] = df["pdb_ids"].fillna("").astype(str).str.strip()
no_pdb = df[df["pdb_ids"] == ""].copy()

OUTDIR.mkdir(parents=True, exist_ok=True)

manifest_rows = []
for _, row in no_pdb.iterrows():
    acc = str(row["uniprot_accession"]).strip()
    orgtype = str(row.get("organism_type", "")).strip()
    gh = str(row.get("gh_family", "")).strip()
    seq = str(row["sequence"]).strip()

    fasta_path = OUTDIR / f"{acc}.fasta"
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(f">{acc}|{orgtype}|{gh}\n{seq}\n")

    manifest_rows.append({
        "uniprot_accession": acc,
        "organism": row.get("organism", ""),
        "organism_type": orgtype,
        "gh_family": gh,
        "length": row.get("length", ""),
        "has_pdb": row.get("has_pdb", 0),
        "input_fasta": str(fasta_path),
        "swiss_model_status": "pending_submission",
        "model_path": "",
        "template": "",
        "gmqe": "",
        "qmean": "",
        "oligomeric_state": "",
        "notes": ""
    })

manifest = pd.DataFrame(manifest_rows)
manifest.to_csv(MANIFEST, index=False)

print(f"Saved manifest: {MANIFEST}")
print(f"Targets without PDB: {len(manifest)}")
