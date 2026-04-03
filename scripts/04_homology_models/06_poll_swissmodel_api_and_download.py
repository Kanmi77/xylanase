#!/usr/bin/env python3
import os
import time
import gzip
import pandas as pd
import requests
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
SUBCSV = BASE / "models" / "swiss_model" / "swissmodel_api_submissions.csv"
OUTCSV = BASE / "models" / "swiss_model" / "swissmodel_api_submissions_polled.csv"
OUTDIR = BASE / "models" / "swiss_model" / "outputs"
LOG = BASE / "results" / "logs" / "swissmodel_api_poll.log"

TOKEN = os.environ.get("SWISSMODEL_TOKEN", "").strip()
if not TOKEN:
    raise SystemExit("SWISSMODEL_TOKEN is not set.")

HEADERS = {"Authorization": f"Token {TOKEN}"}

if OUTCSV.exists():
    df = pd.read_csv(OUTCSV)
else:
    df = pd.read_csv(SUBCSV)

text_cols = [
    "project_id",
    "poll_status",
    "coordinates_url",
    "modelcif_url",
    "downloaded_pdb",
]
for col in text_cols:
    if col not in df.columns:
        df[col] = ""
    df[col] = df[col].fillna("").astype(str)

OUTDIR.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)

def should_poll(status: str, downloaded_pdb: str) -> bool:
    status = str(status).strip()
    downloaded_pdb = str(downloaded_pdb).strip()
    if downloaded_pdb:
        return False
    if status == "FAILED":
        return False
    return True

def fetch_pdb_text(url: str, headers: dict) -> str:
    resp = requests.get(url, headers=headers, timeout=120)
    resp.raise_for_status()

    # API models are returned as .pdb.gz
    if url.endswith(".gz"):
        return gzip.decompress(resp.content).decode("utf-8", errors="replace")

    return resp.text

with open(LOG, "a", encoding="utf-8") as log:
    for idx, r in df.iterrows():
        acc = str(r["uniprot_accession"]).strip()
        project_id = str(r.get("project_id", "")).strip()
        poll_status = str(r.get("poll_status", "")).strip()
        downloaded_pdb = str(r.get("downloaded_pdb", "")).strip()

        if not project_id:
            df.at[idx, "poll_status"] = "NO_PROJECT_ID"
            log.write(f"[SKIP] {acc} no project_id\n")
            continue

        if not should_poll(poll_status, downloaded_pdb):
            continue

        summary_url = f"https://swissmodel.expasy.org/project/{project_id}/models/summary/"
        response_obj = None
        final_status = poll_status if poll_status else ""

        if final_status == "COMPLETED":
            try:
                resp = requests.get(summary_url, headers=HEADERS, timeout=60)
                resp.raise_for_status()
                response_obj = resp.json()
            except Exception as e:
                log.write(f"[ERROR] completed-summary {acc} {e}\n")
        else:
            for attempt in range(12):
                try:
                    resp = requests.get(summary_url, headers=HEADERS, timeout=60)

                    if resp.status_code == 429:
                        wait_s = min(300, 20 * (attempt + 1))
                        final_status = "RATE_LIMITED"
                        log.write(f"[429] {acc} project={project_id} wait={wait_s}s\n")
                        time.sleep(wait_s)
                        continue

                    resp.raise_for_status()
                    response_obj = resp.json()
                    final_status = str(response_obj.get("status", "")).strip()
                    log.write(f"[POLL] {acc} project={project_id} status={final_status}\n")

                    if final_status in ["COMPLETED", "FAILED"]:
                        break

                    time.sleep(15)

                except Exception as e:
                    final_status = f"ERROR: {e}"
                    log.write(f"[ERROR] polling {acc} project={project_id} {e}\n")
                    time.sleep(30)

        if response_obj is not None:
            final_status = str(response_obj.get("status", final_status)).strip()

        df.at[idx, "poll_status"] = final_status

        if final_status == "COMPLETED" and response_obj:
            models = response_obj.get("models", [])
            if models:
                model = models[0]
                coordinates_url = str(model.get("coordinates_url", "")).strip()
                modelcif_url = str(model.get("modelcif_url", "")).strip()

                df.at[idx, "coordinates_url"] = coordinates_url
                df.at[idx, "modelcif_url"] = modelcif_url

                if coordinates_url and not downloaded_pdb:
                    try:
                        pdb_text = fetch_pdb_text(coordinates_url, HEADERS)

                        if pdb_text.startswith("HEADER") or pdb_text.startswith("ATOM") or "ATOM" in pdb_text:
                            out_pdb = OUTDIR / f"{acc}_swissmodel.pdb"
                            out_pdb.write_text(pdb_text, encoding="utf-8")
                            df.at[idx, "downloaded_pdb"] = str(out_pdb)
                            log.write(f"[DOWNLOADED] {acc} -> {out_pdb}\n")
                        else:
                            log.write(f"[BAD_PDB] {acc} {coordinates_url}\n")
                    except Exception as e:
                        log.write(f"[ERROR] download {acc} {e}\n")

        df.to_csv(OUTCSV, index=False)

df.to_csv(OUTCSV, index=False)

print(f"Saved: {OUTCSV}")
print(f"Saved: {LOG}")
print(df["poll_status"].value_counts(dropna=False))
print(f"PDB files in outputs: {len(list(OUTDIR.glob('*.pdb')))}")
