#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import re

BASE = Path.home() / "xylanase" / "xylanase"
MASTER = BASE / "data" / "curated" / "xylanase_master_table_annotated.csv"
FILTERED = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11.csv"
BRENDA = BASE / "data" / "brenda" / "processed" / "brenda_xylanase_structured.csv"

MASTER_OUT = BASE / "data" / "curated" / "xylanase_master_table_with_brenda.csv"
FILTERED_OUT = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv"
REPORT = BASE / "results" / "reports" / "brenda_merge_report.txt"

def norm_text(x):
    if pd.isna(x):
        return ""
    x = str(x).strip().lower()
    x = re.sub(r"\s+", " ", x)
    return x

def make_wide(brenda):
    wide_parts = []

    for record_type, prefix in [
        ("optimum_temperature", "optimum_temperature"),
        ("temperature_range", "temperature_range"),
        ("temperature_stability", "temperature_stability"),
        ("optimum_pH", "optimum_pH"),
        ("pH_range", "pH_range"),
    ]:
        sub = brenda[brenda["record_type"] == record_type].copy()
        if sub.empty:
            continue

        sub = sub.rename(columns={
            "value": prefix,
            "commentary": f"{prefix}_commentary",
            "literature": f"{prefix}_literature",
            "source_file": f"{prefix}_source_file",
        })

        keep = [
            "uniprot_accession",
            "organism_norm",
            prefix,
            f"{prefix}_commentary",
            f"{prefix}_literature",
            f"{prefix}_source_file",
        ]
        sub = sub[keep].drop_duplicates(subset=["uniprot_accession", "organism_norm"])
        wide_parts.append(sub)

    if not wide_parts:
        return pd.DataFrame(columns=["uniprot_accession", "organism_norm"])

    wide = wide_parts[0]
    for part in wide_parts[1:]:
        wide = wide.merge(part, how="outer", on=["uniprot_accession", "organism_norm"])

    return wide

def merge_dataset(df, wide):
    df = df.copy()
    df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()
    df["organism_norm"] = df["organism"].map(norm_text)

    wide = wide.copy()
    wide["uniprot_accession"] = wide["uniprot_accession"].fillna("").astype(str).str.strip()
    wide["organism_norm"] = wide["organism_norm"].fillna("").astype(str)

    up_wide = wide[wide["uniprot_accession"] != ""].drop_duplicates(subset=["uniprot_accession"])
    out = df.merge(up_wide.drop(columns=["organism_norm"]), how="left", on="uniprot_accession")

    org_wide = wide.drop_duplicates(subset=["organism_norm"])
    out = out.merge(org_wide, how="left", on="organism_norm", suffixes=("", "_org"))

    for col in [c for c in out.columns if c.endswith("_org")]:
        base = col[:-4]
        out[base] = out[base].combine_first(out[col])
        out = out.drop(columns=[col])

    if "optimum_temperature" in out.columns and "optimum_temperature_brenda" in out.columns:
        out["optimum_temperature"] = out["optimum_temperature"].fillna("").astype(str)
        out["optimum_temperature_brenda"] = out["optimum_temperature_brenda"].fillna("").astype(str)
        out["optimum_temperature"] = out.apply(
            lambda r: r["optimum_temperature"] if r["optimum_temperature"].strip() else r["optimum_temperature_brenda"],
            axis=1
        )

    if "optimum_pH" in out.columns and "optimum_pH_brenda" in out.columns:
        out["optimum_pH"] = out["optimum_pH"].fillna("").astype(str)
        out["optimum_pH_brenda"] = out["optimum_pH_brenda"].fillna("").astype(str)
        out["optimum_pH"] = out.apply(
            lambda r: r["optimum_pH"] if r["optimum_pH"].strip() else r["optimum_pH_brenda"],
            axis=1
        )

    out = out.drop(columns=["organism_norm"])
    return out

master = pd.read_csv(MASTER)
filtered = pd.read_csv(FILTERED)
brenda = pd.read_csv(BRENDA)

brenda["uniprot_accession"] = brenda["uniprot_accession"].fillna("").astype(str).str.strip()
brenda["organism_norm"] = brenda["organism"].map(norm_text)

wide = make_wide(brenda)

rename_fill = {}
if "optimum_temperature" in wide.columns:
    rename_fill["optimum_temperature"] = "optimum_temperature_brenda"
if "optimum_pH" in wide.columns:
    rename_fill["optimum_pH"] = "optimum_pH_brenda"
wide = wide.rename(columns=rename_fill)

master_out = merge_dataset(master, wide)
filtered_out = merge_dataset(filtered, wide)

master_out.to_csv(MASTER_OUT, index=False)
filtered_out.to_csv(FILTERED_OUT, index=False)

master_temp = (master_out["optimum_temperature"].fillna("").astype(str).str.strip() != "").sum()
filtered_temp = (filtered_out["optimum_temperature"].fillna("").astype(str).str.strip() != "").sum()

master_ph = (master_out["optimum_pH"].fillna("").astype(str).str.strip() != "").sum()
filtered_ph = (filtered_out["optimum_pH"].fillna("").astype(str).str.strip() != "").sum()

lines = [
    "BRENDA merge report",
    "=" * 50,
    f"Master rows: {len(master_out)}",
    f"Filtered rows: {len(filtered_out)}",
    f"Master rows with optimum_temperature after merge: {master_temp}",
    f"Filtered rows with optimum_temperature after merge: {filtered_temp}",
    f"Master rows with optimum_pH after merge: {master_ph}",
    f"Filtered rows with optimum_pH after merge: {filtered_ph}",
]

REPORT.write_text("\n".join(lines), encoding="utf-8")

print("\n".join(lines))
print(f"Saved: {MASTER_OUT}")
print(f"Saved: {FILTERED_OUT}")
print(f"Saved: {REPORT}")
