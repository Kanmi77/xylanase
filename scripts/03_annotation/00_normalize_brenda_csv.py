#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "brenda" / "raw" / "brenda_xylanase_optima.csv"
OUTFILE = BASE / "data" / "brenda" / "processed" / "brenda_xylanase_structured.csv"

df = pd.read_csv(INFILE)

# expected flexible columns from your manual export or cleanup
rename_map = {
    "type": "record_type",
    "record_type": "record_type",
    "organism": "organism",
    "uniprot": "uniprot_accession",
    "uniprot_accession": "uniprot_accession",
    "value": "value",
    "temperature": "value",
    "commentary": "commentary",
    "comment": "commentary",
    "literature": "literature",
    "source_pdf": "source_pdf"
}

df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

required = ["record_type", "organism", "uniprot_accession", "value", "commentary", "literature", "source_pdf"]
for col in required:
    if col not in df.columns:
        df[col] = ""

df = df[required].copy()
df["record_type"] = df["record_type"].astype(str).str.strip()
df["organism"] = df["organism"].astype(str).str.strip()
df["uniprot_accession"] = df["uniprot_accession"].astype(str).str.strip()
df["value"] = df["value"].astype(str).str.strip()
df["commentary"] = df["commentary"].astype(str).str.strip()
df["literature"] = df["literature"].astype(str).str.strip()
df["source_pdf"] = df["source_pdf"].astype(str).str.strip()

OUTFILE.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(df)}")
print(df["record_type"].value_counts(dropna=False))
