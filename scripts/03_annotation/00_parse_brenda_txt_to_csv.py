#!/usr/bin/env python3
import re
import csv
from pathlib import Path

BASE = Path("/home/ubuntu/xylanase/xylanase/data/brenda/raw")

FILES = [
    ("brenda.txt",  "temperature_optimum",  "brenda_temperature_optimum.csv"),
    ("brenda2.txt", "temperature_range",    "brenda_temperature_range.csv"),
    ("brenda3.txt", "temperature_stability","brenda_temperature_stability.csv"),
]

def clean_lines(path):
    lines = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("2/") or s.startswith("about:blank"):
            continue
        if s.startswith("Information on EC 3.2.1.8"):
            continue
        if s in {"TEMPERATURE ORGANISM", "UNIPROT", "COMMENTARY", "LITERATURE"}:
            continue
        if s in {"OPTIMUM", "RANGE", "STABILITY"}:
            continue
        lines.append(s)
    return lines

def is_temp_or_range(s):
    return bool(re.fullmatch(r"\d+(?:\s*-\s*\d+)?", s))

def is_uniprotish(s):
    if s == "-":
        return True
    parts = [x.strip() for x in s.split(",")]
    for p in parts:
        if not p:
            continue
        if re.fullmatch(r"[A-Z0-9]{6,10}", p) is None:
            return False
    return True

def is_literature(s):
    return bool(re.fullmatch(r"\d{5,7}", s))

def parse_records(lines, section):
    records = []
    i = 0
    n = len(lines)

    while i < n:
        if not is_temp_or_range(lines[i]):
            i += 1
            continue

        temp_value = lines[i]
        i += 1

        organism = []
        while i < n and not is_uniprotish(lines[i]):
            if is_literature(lines[i]):
                break
            organism.append(lines[i])
            i += 1

        uniprot = ""
        if i < n and is_uniprotish(lines[i]):
            uniprot = lines[i]
            i += 1

        commentary = []
        while i < n and not is_literature(lines[i]):
            if is_temp_or_range(lines[i]):
                break
            commentary.append(lines[i])
            i += 1

        literature = ""
        if i < n and is_literature(lines[i]):
            literature = lines[i]
            i += 1

        records.append({
            "section": section,
            "temperature_value": temp_value,
            "organism": " ".join(organism).strip(),
            "uniprot": uniprot.strip(),
            "commentary": " ".join(commentary).strip(),
            "literature": literature.strip(),
        })

    return records

for infile, section, outfile in FILES:
    path = BASE / infile
    lines = clean_lines(path)
    records = parse_records(lines, section)

    outpath = BASE.parent / "processed" / outfile
    outpath.parent.mkdir(parents=True, exist_ok=True)

    with open(outpath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["section", "temperature_value", "organism", "uniprot", "commentary", "literature"]
        )
        writer.writeheader()
        writer.writerows(records)

    print(f"Saved: {outpath} ({len(records)} rows)")
