#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"

FILES = [
    BASE / "data" / "curated" / "xylanase_master_table_with_brenda.csv",
    BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11_with_brenda.csv",
]

DROP_PATTERNS = [
    "_brenda",
    "_commentary",
    "_literature",
    "_source_file",
]

DROP_EXACT = {
    "temperature_range",
    "pH_range",
}

for f in FILES:
    if not f.exists():
        print(f"[MISSING] {f}")
        continue

    df = pd.read_csv(f)
    original_cols = df.columns.tolist()

    keep_cols = []
    removed_cols = []

    for c in original_cols:
        if c in DROP_EXACT or any(p in c for p in DROP_PATTERNS):
            removed_cols.append(c)
        else:
            keep_cols.append(c)

    cleaned = df[keep_cols].copy()
    cleaned.to_csv(f, index=False)

    print(f"\nCleaned: {f}")
    print(f"Original columns: {len(original_cols)}")
    print(f"Kept columns: {len(keep_cols)}")
    print("Removed:")
    print(removed_cols)
