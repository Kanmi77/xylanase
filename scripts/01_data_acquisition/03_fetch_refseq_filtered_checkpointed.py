#!/usr/bin/env python3
import time
import requests
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "data" / "curated" / "xylanase_filtered_bacterial_fungal_GH10_GH11.csv"
OUT_FASTA = BASE / "data" / "refseq" / "processed" / "filtered_refseq_sequences.fasta"
OUT_META = BASE / "data" / "refseq" / "processed" / "filtered_refseq_metadata.csv"
OUT_DONE = BASE / "data" / "refseq" / "processed" / "filtered_refseq_done.txt"
OUT_FAIL = BASE / "data" / "refseq" / "processed" / "filtered_refseq_failed.txt"

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EMAIL = "sunkanmikamal.edu@gmail.com"
TOOL = "xylanase_filtered_refseq"

def split_ids(value):
    if pd.isna(value):
        return []
    return [x.strip().split()[0] for x in str(value).replace(";", ",").split(",") if x.strip()]

df = pd.read_csv(INFILE)

unique_ids = set()
for v in df["refseq_ids"]:
    unique_ids.update(split_ids(v))

unique_ids = sorted(unique_ids)

done = set()
if OUT_DONE.exists():
    done = set(x.strip() for x in OUT_DONE.read_text().splitlines() if x.strip())

failed = set()
if OUT_FAIL.exists():
    failed = set(x.strip() for x in OUT_FAIL.read_text().splitlines() if x.strip())

remaining = [x for x in unique_ids if x not in done and x not in failed]

print(f"Total unique RefSeq IDs in filtered subset: {len(unique_ids)}")
print(f"Already done: {len(done)}")
print(f"Already failed: {len(failed)}")
print(f"Remaining: {len(remaining)}")

OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)

with open(OUT_FASTA, "a", encoding="utf-8") as fasta_out, \
     open(OUT_DONE, "a", encoding="utf-8") as done_out, \
     open(OUT_FAIL, "a", encoding="utf-8") as fail_out:

    for i, rid in enumerate(remaining, start=1):
        try:
            r = requests.get(
                EFETCH,
                params={
                    "db": "protein",
                    "id": rid,
                    "rettype": "fasta",
                    "retmode": "text",
                    "tool": TOOL,
                    "email": EMAIL,
                },
                timeout=60,
            )
            r.raise_for_status()
            txt = r.text.strip()
            if txt.startswith(">"):
                fasta_out.write(txt + "\n")
                done_out.write(rid + "\n")
                print(f"[OK] {i}/{len(remaining)} {rid}")
            else:
                fail_out.write(rid + "\n")
                print(f"[NO FASTA] {i}/{len(remaining)} {rid}")
        except Exception:
            fail_out.write(rid + "\n")
            print(f"[FAILED] {i}/{len(remaining)} {rid}")

        fasta_out.flush()
        done_out.flush()
        fail_out.flush()
        time.sleep(0.34)
