#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
RAW = BASE / "data" / "brenda" / "raw"
OUT = BASE / "data" / "brenda" / "processed" / "brenda_xylanase_structured.csv"

FILES = [
    ("xylanase_temperature_optimum_FIXED.csv", "optimum_temperature"),
    ("xylanase_temperature_range_FIXED.csv", "temperature_range"),
    ("xylanase_temperature_stability_FIXED.csv", "temperature_stability"),
    ("xylanase_ph_optimum.csv", "optimum_pH"),
    ("xylanase_ph_range.csv", "pH_range"),
]

STANDARD_COLS = [
    "record_type",
    "value",
    "organism",
    "uniprot_accession",
    "commentary",
    "literature",
    "source_file",
]

def normalize_colname(col: str) -> str:
    return str(col).strip().lower().replace(" ", "_")

def standardize_columns(df: pd.DataFrame, record_type: str, source_file: str) -> pd.DataFrame:
    original_cols = df.columns.tolist()
    df = df.copy()
    df.columns = [normalize_colname(c) for c in df.columns]

    rename_map = {
        "temperature_value": "value",
        "ph_value": "value",
	"phoptimum": "value",
        "ph_optimum": "value",
        "ph_range": "value",
        "range_value": "value",
        "value": "value",
        "temperature": "value",
        "ph": "value",
        "uniprot": "uniprot_accession",
        "uniprot_accession": "uniprot_accession",
        "organism": "organism",
        "comment": "commentary",
        "comments": "commentary",
        "commentary": "commentary",
        "literature": "literature",
        "reference": "literature",
        "references": "literature",
    }

    matched = {}
    for c in df.columns:
        if c in rename_map:
            matched[c] = rename_map[c]
    df = df.rename(columns=matched)

    for col in ["value", "organism", "uniprot_accession", "commentary", "literature"]:
        if col not in df.columns:
            df[col] = ""

    out = df[["value", "organism", "uniprot_accession", "commentary", "literature"]].copy()
    for c in out.columns:
        out[c] = out[c].fillna("").astype(str).str.strip()

    out["record_type"] = record_type
    out["source_file"] = source_file
    out = out[STANDARD_COLS]

    print(f"\n[{source_file}]")
    print("Original columns:", original_cols)
    print("Normalized columns:", df.columns.tolist())
    print("Non-empty counts:")
    for c in ["value", "organism", "uniprot_accession", "commentary", "literature"]:
        print(f"  {c}: {(out[c] != '').sum()}")

    return out

all_frames = []

for filename, record_type in FILES:
    path = RAW / filename
    if not path.exists():
        print(f"[MISSING] {path}")
        continue

    df = pd.read_csv(path)
    clean = standardize_columns(df, record_type, filename)
    all_frames.append(clean)
    print(f"[OK] {filename}: {len(clean)} rows")

if not all_frames:
    raise SystemExit("No source BRENDA files found.")

merged = pd.concat(all_frames, ignore_index=True)

# remove only truly empty rows
merged = merged[
    ~(
        (merged["value"] == "") &
        (merged["organism"] == "") &
        (merged["uniprot_accession"] == "") &
        (merged["commentary"] == "") &
        (merged["literature"] == "")
    )
].copy()

# deduplicate only on meaningful content
merged = merged.drop_duplicates(
    subset=["record_type", "value", "organism", "uniprot_accession", "commentary", "literature", "source_file"]
)

OUT.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(OUT, index=False)

print(f"\nSaved: {OUT}")
print(f"Total rows: {len(merged)}")
print("\nBy record type:")
print(merged["record_type"].value_counts(dropna=False))
