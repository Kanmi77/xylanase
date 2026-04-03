#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
if not INFILE.exists():
    INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11.csv"

OUTDIR = BASE / "models" / "swiss_model" / "inputs"
MANIFEST = BASE / "models" / "swiss_model" / "swiss_model_manifest_final.csv"

df = pd.read_csv(INFILE)
df["pdb_ids"] = df["pdb_ids"].fillna("").astype(str).str.strip()

no_pdb = df[df["pdb_ids"] == ""].copy()

OUTDIR.mkdir(parents=True, exist_ok=True)

rows = []
for _, r in no_pdb.iterrows():
    acc = str(r["uniprot_accession"]).strip()
    seq = str(r["sequence"]).strip()
    org = str(r.get("organism", "")).strip()
    org_type = str(r.get("organism_type", "")).strip()
    gh = str(r.get("gh_family", "")).strip()

    fasta_path = OUTDIR / f"{acc}.fasta"
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(f">{acc}|{org_type}|{gh}\n{seq}\n")

    rows.append({
        "uniprot_accession": acc,
        "organism": org,
        "organism_type": org_type,
        "gh_family": gh,
        "length": r.get("length", ""),
        "input_fasta": str(fasta_path),
        "model_source": "",
        "swiss_model_status": "pending_repository_check",
        "model_path": "",
        "template": "",
        "gmqe": "",
        "qmean": "",
        "notes": ""
    })

manifest = pd.DataFrame(rows)
manifest.to_csv(MANIFEST, index=False)

print(f"Saved: {MANIFEST}")
print(f"Targets without PDB: {len(manifest)}")
