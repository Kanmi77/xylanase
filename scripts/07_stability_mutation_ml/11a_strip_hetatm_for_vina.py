#!/usr/bin/env python3
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INDIR = BASE / "models" / "swiss_model" / "outputs"
OUTDIR = BASE / "results" / "docking" / "cleaned_receptors"

OUTDIR.mkdir(parents=True, exist_ok=True)

count = 0
for pdb in INDIR.glob("*.pdb"):
    out = OUTDIR / pdb.name
    kept = []
    with open(pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec in {"ATOM", "TER", "END", "ENDMDL", "MODEL"}:
                kept.append(line)
    with open(out, "w", encoding="utf-8") as f:
        f.writelines(kept)
    count += 1

print(f"Saved cleaned PDBs to: {OUTDIR}")
print(f"Count: {count}")
