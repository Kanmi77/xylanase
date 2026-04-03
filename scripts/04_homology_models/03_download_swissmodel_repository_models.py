#!/usr/bin/env python3
import time
import requests
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
MANIFEST = BASE / "models" / "swiss_model" / "swiss_model_manifest_final.csv"
OUTDIR = BASE / "models" / "swiss_model" / "outputs"
OUTCSV = BASE / "models" / "swiss_model" / "swiss_model_repository_hits.csv"
LOG = BASE / "results" / "logs" / "swiss_model_repository_download.log"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)

HEADERS = {
    "User-Agent": "xylanase-thesis-swissmodel-check/1.0"
}

def get_best_structure(entry_json):
    result = entry_json.get("result", {})
    structures = result.get("structures", [])
    if not structures:
        return None
    # choose best by highest GMQE, then highest identity
    def score(s):
        gmqe = s.get("gmqe", -1)
        identity = s.get("identity", -1)
        coverage = s.get("coverage", -1)
        return (gmqe, identity, coverage)
    return sorted(structures, key=score, reverse=True)[0]

def main():
    manifest = pd.read_csv(MANIFEST)
    rows = []

    with open(LOG, "w", encoding="utf-8") as log:
        for _, r in manifest.iterrows():
            acc = str(r["uniprot_accession"]).strip()
            if not acc:
                continue

            out_pdb = OUTDIR / f"{acc}_swissmodel.pdb"
            api_url = f"https://swissmodel.expasy.org/repository/uniprot/{acc}.json?provider=swissmodel"

            try:
                resp = requests.get(api_url, headers=HEADERS, timeout=60)
                resp.raise_for_status()
                data = resp.json()

                best = get_best_structure(data)
                if best is None:
                    log.write(f"[NO_MODEL] {acc}\n")
                    rows.append({
                        "uniprot_accession": acc,
                        "repository_status": "no_model",
                        "coordinates_url": "",
                        "template": "",
                        "gmqe": "",
                        "identity": "",
                        "coverage": "",
                        "qmean4_z_score": "",
                        "downloaded_pdb": ""
                    })
                    time.sleep(0.2)
                    continue

                pdb_url = best.get("coordinates", "")
                template = best.get("template", "")
                gmqe = best.get("gmqe", "")
                identity = best.get("identity", "")
                coverage = best.get("coverage", "")
                qmean4_z = best.get("qmean", {}).get("qmean4_z_score", "")

                downloaded = ""
                if pdb_url:
                    pdb_resp = requests.get(pdb_url, headers=HEADERS, timeout=60)
                    pdb_resp.raise_for_status()
                    text = pdb_resp.text
                    if text.startswith("HEADER") or text.startswith("ATOM") or "ATOM" in text:
                        out_pdb.write_text(text, encoding="utf-8")
                        downloaded = str(out_pdb)
                        log.write(f"[DOWNLOADED] {acc} -> {out_pdb.name}\n")
                    else:
                        log.write(f"[BAD_PDB] {acc} url={pdb_url}\n")
                else:
                    log.write(f"[NO_COORD_URL] {acc}\n")

                rows.append({
                    "uniprot_accession": acc,
                    "repository_status": "model_found" if pdb_url else "model_found_no_coordinates",
                    "coordinates_url": pdb_url,
                    "template": template,
                    "gmqe": gmqe,
                    "identity": identity,
                    "coverage": coverage,
                    "qmean4_z_score": qmean4_z,
                    "downloaded_pdb": downloaded
                })

            except Exception as e:
                log.write(f"[ERROR] {acc} {e}\n")
                rows.append({
                    "uniprot_accession": acc,
                    "repository_status": f"error: {e}",
                    "coordinates_url": "",
                    "template": "",
                    "gmqe": "",
                    "identity": "",
                    "coverage": "",
                    "qmean4_z_score": "",
                    "downloaded_pdb": ""
                })

            time.sleep(0.2)

    out = pd.DataFrame(rows)
    out.to_csv(OUTCSV, index=False)

    print(f"Saved: {OUTCSV}")
    print(f"Saved: {LOG}")
    print(out['repository_status'].value_counts(dropna=False))
    print(f"PDB files in outputs: {len(list(OUTDIR.glob('*.pdb')))}")

if __name__ == "__main__":
    main()
