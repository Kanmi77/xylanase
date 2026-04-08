#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
OUTFILE = BASE / "results" / "features" / "protparam_features.csv"

df = pd.read_csv(INFILE).copy()
df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()
df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()

aa_order = list("ACDEFGHIKLMNPQRSTVWY")
rows = []

for _, row in df.iterrows():
    acc = row["uniprot_accession"]
    seq = row["sequence"]

    if not acc or not seq:
        continue

    try:
        pa = ProteinAnalysis(seq)
        aa_percent = pa.amino_acids_percent
        ss = pa.secondary_structure_fraction()

        rec = {
            "uniprot_accession": acc,
            "organism": row.get("organism", ""),
            "organism_type": row.get("organism_type", ""),
            "gh_family": row.get("gh_family", ""),
            "length": len(seq),
            "molecular_weight": pa.molecular_weight(),
            "aromaticity": pa.aromaticity(),
            "instability_index": pa.instability_index(),
            "gravy": pa.gravy(),
            "predicted_pI": pa.isoelectric_point(),
            "secondary_structure_fraction_helix": ss[0],
            "secondary_structure_fraction_turn": ss[1],
            "secondary_structure_fraction_sheet": ss[2],
            "optimum_temperature": row.get("optimum_temperature", ""),
            "optimum_pH": row.get("optimum_pH", ""),
        }

        for aa in aa_order:
            rec[f"aa_{aa}_fraction"] = aa_percent.get(aa, 0.0)

        rows.append(rec)

    except Exception as e:
        rows.append({
            "uniprot_accession": acc,
            "organism": row.get("organism", ""),
            "organism_type": row.get("organism_type", ""),
            "gh_family": row.get("gh_family", ""),
            "length": len(seq),
            "protparam_error": str(e),
        })

out = pd.DataFrame(rows)
OUTFILE.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(out)}")
