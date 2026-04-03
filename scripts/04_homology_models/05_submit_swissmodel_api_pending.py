#!/usr/bin/env python3
import os
import time
import pandas as pd
import requests
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
PENDING = BASE / "models" / "swiss_model" / "pending_submission_targets.csv"
OUTCSV = BASE / "models" / "swiss_model" / "swissmodel_api_submissions.csv"
LOG = BASE / "results" / "logs" / "swissmodel_api_submit.log"

TOKEN = os.environ.get("SWISSMODEL_TOKEN", "").strip()
if not TOKEN:
    raise SystemExit("SWISSMODEL_TOKEN is not set.")

HEADERS = {"Authorization": f"Token {TOKEN}"}
URL = "https://swissmodel.expasy.org/automodel"

df = pd.read_csv(PENDING)
rows = []

LOG.parent.mkdir(parents=True, exist_ok=True)

with open(LOG, "w", encoding="utf-8") as log:
    for i, r in df.iterrows():
        acc = str(r["uniprot_accession"]).strip()
        fasta_path = str(r["input_fasta"]).strip()

        if not acc or not fasta_path:
            log.write(f"[SKIP] missing accession or fasta at row {i}\n")
            continue

        seq = []
        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.startswith(">"):
                    seq.append(line.strip())
        sequence = "".join(seq)

        payload = {
            "target_sequences": [sequence],
            "project_title": f"xylanase_{acc}"
        }

        try:
            resp = requests.post(URL, headers=HEADERS, json=payload, timeout=120)
            status_code = resp.status_code
            ok = resp.ok

            try:
                data = resp.json()
            except Exception:
                data = {"raw_text": resp.text[:1000]}

            project_id = data.get("project_id", "")
            rows.append({
                "uniprot_accession": acc,
                "project_title": f"xylanase_{acc}",
                "project_id": project_id,
                "submit_status_code": status_code,
                "submit_ok": ok,
                "submit_response": str(data)[:5000],
                "poll_status": "",
                "coordinates_url": "",
                "modelcif_url": "",
                "downloaded_pdb": ""
            })

            log.write(f"[SUBMIT] {acc} status={status_code} project_id={project_id}\n")

        except Exception as e:
            rows.append({
                "uniprot_accession": acc,
                "project_title": f"xylanase_{acc}",
                "project_id": "",
                "submit_status_code": "",
                "submit_ok": False,
                "submit_response": str(e),
                "poll_status": "",
                "coordinates_url": "",
                "modelcif_url": "",
                "downloaded_pdb": ""
            })
            log.write(f"[ERROR] {acc} {e}\n")

        time.sleep(0.5)

out = pd.DataFrame(rows)
out.to_csv(OUTCSV, index=False)

print(f"Saved: {OUTCSV}")
print(f"Saved: {LOG}")
print(out[['submit_status_code','submit_ok']].value_counts(dropna=False))
