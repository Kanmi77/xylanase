#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import re

BASE = Path.home() / "xylanase" / "xylanase"

MASTER = BASE / "data" / "curated" / "xylanase_master_table_annotated.csv"
if not MASTER.exists():
    MASTER = BASE / "data" / "curated" / "xylanase_master_table.csv"

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

def prepare_brenda_wide(brenda):
    opt = brenda[brenda["record_type"] == "optimum_temperature"].copy()
    opt["optimum_temperature_brenda"] = opt["temperature_value"]
    opt["optimum_temperature_source"] = "BRENDA:" + opt["literature"].fillna("").astype(str).str.strip()
    rng = brenda[brenda["record_type"] == "temperature_range"].copy()
    rng["temperature_range_brenda"] = rng["temperature_value"]
    rng["temperature_range_source"] = "BRENDA:" + rng["literature"].fillna("").astype(str).str.strip()
    stab = brenda[brenda["record_type"] == "temperature_stability"].copy()
    stab["temperature_stability_brenda"] = stab["temperature_value"]
    stab["temperature_stability_source"] = "BRENDA:" + stab["literature"].fillna("").astype(str).str.strip()

    opt = opt.rename(columns={
        "value": "optimum_temperature_brenda",
        "commentary": "optimum_temperature_commentary",
        "literature": "optimum_temperature_literature",
        "source_pdf": "optimum_temperature_source"
    })

    rng = rng.rename(columns={
        "value": "temperature_range_brenda",
        "commentary": "temperature_range_commentary",
        "literature": "temperature_range_literature",
        "source_pdf": "temperature_range_source"
    })

    stab = stab.rename(columns={
        "value": "temperature_stability_brenda",
        "commentary": "temperature_stability_commentary",
        "literature": "temperature_stability_literature",
        "source_pdf": "temperature_stability_source"
    })

    keep_opt = ["uniprot_accession", "organism_norm", "optimum_temperature_brenda",
                "optimum_temperature_commentary", "optimum_temperature_literature",
                "optimum_temperature_source"]
    keep_rng = ["uniprot_accession", "organism_norm", "temperature_range_brenda",
                "temperature_range_commentary", "temperature_range_literature",
                "temperature_range_source"]
    keep_stab = ["uniprot_accession", "organism_norm", "temperature_stability_brenda",
                 "temperature_stability_commentary", "temperature_stability_literature",
                 "temperature_stability_source"]

    opt = opt[keep_opt].drop_duplicates(subset=["uniprot_accession", "organism_norm"])
    rng = rng[keep_rng].drop_duplicates(subset=["uniprot_accession", "organism_norm"])
    stab = stab[keep_stab].drop_duplicates(subset=["uniprot_accession", "organism_norm"])

    wide = opt.merge(rng, how="outer", on=["uniprot_accession", "organism_norm"])
    wide = wide.merge(stab, how="outer", on=["uniprot_accession", "organism_norm"])
    return wide

def merge_dataset(df, wide):
    df = df.copy()
    df["uniprot_accession"] = df["uniprot_accession"].fillna("").astype(str).str.strip()
    df["organism_norm"] = df["organism"].map(norm_text)

    wide = wide.copy()
    wide["uniprot_accession"] = wide["uniprot_accession"].fillna("").astype(str).str.strip()
    wide["organism_norm"] = wide["organism_norm"].fillna("").astype(str)

    # first exact UniProt match
    up_wide = wide[wide["uniprot_accession"] != ""].drop_duplicates(subset=["uniprot_accession"])
    out = df.merge(
        up_wide.drop(columns=["organism_norm"]),
        how="left",
        on="uniprot_accession"
    )

    # second organism fallback only for rows still missing optimum temperature
    org_wide = wide.drop_duplicates(subset=["organism_norm"])
    out = out.merge(
        org_wide,
        how="left",
        on="organism_norm",
        suffixes=("", "_org")
    )

    # fill optimum temperature from organism fallback only if exact UniProt gave nothing
    for col in [
        "optimum_temperature_brenda",
        "optimum_temperature_commentary",
        "optimum_temperature_literature",
        "optimum_temperature_source",
        "temperature_range_brenda",
        "temperature_range_commentary",
        "temperature_range_literature",
        "temperature_range_source",
        "temperature_stability_brenda",
        "temperature_stability_commentary",
        "temperature_stability_literature",
        "temperature_stability_source",
    ]:
        org_col = f"{col}_org"
        if org_col in out.columns:
            out[col] = out[col].combine_first(out[org_col])
            out = out.drop(columns=[org_col])

    # map optimum temperature into main field only if currently empty
    if "optimum_temperature" in out.columns:
        out["optimum_temperature"] = out["optimum_temperature"].fillna("")
        out["optimum_temperature_brenda"] = out["optimum_temperature_brenda"].fillna("")
        out["optimum_temperature"] = out.apply(
            lambda r: r["optimum_temperature"]
            if str(r["optimum_temperature"]).strip() != ""
            else r["optimum_temperature_brenda"],
            axis=1
        )

    out = out.drop(columns=["organism_norm"])
    return out

master = pd.read_csv(MASTER)
filtered = pd.read_csv(FILTERED)
brenda = pd.read_csv(BRENDA)

brenda["uniprot_accession"] = brenda["uniprot_accession"].fillna("").astype(str).str.strip()
brenda["organism_norm"] = brenda["organism"].map(norm_text)

wide = prepare_brenda_wide(brenda)

master_out = merge_dataset(master, wide)
filtered_out = merge_dataset(filtered, wide)

master_out.to_csv(MASTER_OUT, index=False)
filtered_out.to_csv(FILTERED_OUT, index=False)

master_temp = (master_out["optimum_temperature"].fillna("").astype(str).str.strip() != "").sum()
filtered_temp = (filtered_out["optimum_temperature"].fillna("").astype(str).str.strip() != "").sum()

lines = []
lines.append("BRENDA merge report")
lines.append("=" * 50)
lines.append(f"Master rows: {len(master_out)}")
lines.append(f"Filtered rows: {len(filtered_out)}")
lines.append(f"Master rows with optimum_temperature after merge: {master_temp}")
lines.append(f"Filtered rows with optimum_temperature after merge: {filtered_temp}")

REPORT.write_text("\n".join(lines), encoding="utf-8")

print("\n".join(lines))
print(f"Saved: {MASTER_OUT}")
print(f"Saved: {FILTERED_OUT}")
print(f"Saved: {REPORT}")
