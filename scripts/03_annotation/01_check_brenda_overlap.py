#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import re

BASE = Path.home() / "xylanase" / "xylanase"
MASTER = BASE / "data" / "curated" / "xylanase_master_table_annotated.csv"
FILTERED = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11.csv"
BRENDA = BASE / "data" / "brenda" / "processed" / "brenda_xylanase_structured.csv"
OUT = BASE / "results" / "reports" / "brenda_overlap_report.txt"

def norm_text(x):
    if pd.isna(x):
        return ""
    x = str(x).strip().lower()
    x = re.sub(r"\s+", " ", x)
    return x

master = pd.read_csv(MASTER)
filtered = pd.read_csv(FILTERED)
brenda = pd.read_csv(BRENDA)

for df in [master, filtered, brenda]:
    if "uniprot_accession" not in df.columns:
        df["uniprot_accession"] = ""
    if "organism" not in df.columns:
        df["organism"] = ""

master["uniprot_accession"] = master["uniprot_accession"].fillna("").astype(str).str.strip()
filtered["uniprot_accession"] = filtered["uniprot_accession"].fillna("").astype(str).str.strip()
brenda["uniprot_accession"] = brenda["uniprot_accession"].fillna("").astype(str).str.strip()

master["organism_norm"] = master["organism"].map(norm_text)
filtered["organism_norm"] = filtered["organism"].map(norm_text)
brenda["organism_norm"] = brenda["organism"].map(norm_text)

master_up = set(master["uniprot_accession"]) - {""}
filtered_up = set(filtered["uniprot_accession"]) - {""}
brenda_up = set(brenda["uniprot_accession"]) - {""}

master_org = set(master["organism_norm"]) - {""}
filtered_org = set(filtered["organism_norm"]) - {""}
brenda_org = set(brenda["organism_norm"]) - {""}

master_up_overlap = sorted(master_up & brenda_up)
filtered_up_overlap = sorted(filtered_up & brenda_up)
master_org_overlap = sorted(master_org & brenda_org)
filtered_org_overlap = sorted(filtered_org & brenda_org)

lines = [
    "BRENDA overlap report",
    "=" * 50,
    f"Master rows: {len(master)}",
    f"Filtered rows: {len(filtered)}",
    f"BRENDA rows: {len(brenda)}",
    "",
    f"Master UniProt overlap: {len(master_up_overlap)}",
    f"Filtered UniProt overlap: {len(filtered_up_overlap)}",
    f"Master organism overlap: {len(master_org_overlap)}",
    f"Filtered organism overlap: {len(filtered_org_overlap)}",
    "",
    "Sample filtered UniProt overlaps:",
]
lines.extend(filtered_up_overlap[:30] if filtered_up_overlap else ["None"])
lines.append("")
lines.append("Sample filtered organism overlaps:")
lines.extend(filtered_org_overlap[:30] if filtered_org_overlap else ["None"])

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(lines), encoding="utf-8")

print("\n".join(lines))
print(f"\nSaved: {OUT}")
