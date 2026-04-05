#!/usr/bin/env python3
import docx
import pandas as pd
import re
from pathlib import Path

DOCX_PATH = Path("results/mutations/Duet_results.docx")
OUT_PATH = Path("results/mutations/duet_results_manual.csv")

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

HEADER_RE = re.compile(r"^([A-Za-z0-9_]+)\s*-\s*([A-Z]\d+[A-Z])\s*$")

def clean_text(text: str) -> str:
    text = text.replace("\xa0", " ")
    text = text.replace("–", "-").replace("—", "-")
    text = text.replace("`", "")
    text = re.sub(r"[ \t]+", " ", text)
    return text.strip()

def load_docx_text(path: Path) -> str:
    doc = docx.Document(str(path))
    paras = [clean_text(p.text) for p in doc.paragraphs]
    paras = [p for p in paras if p]
    return "\n".join(paras)

def split_into_blocks(text: str) -> list[str]:
    lines = [clean_text(x) for x in text.splitlines()]
    lines = [x for x in lines if x]

    blocks = []
    current = []

    for line in lines:
        if HEADER_RE.match(line):
            if current:
                blocks.append("\n".join(current))
            current = [line]
        else:
            if current:
                current.append(line)

    if current:
        blocks.append("\n".join(current))

    return blocks

def extract_first_float(text: str):
    m = re.search(r"([+-]?\d+(?:\.\d+)?)", text)
    return float(m.group(1)) if m else None

def extract_ddg(block: str, tool_name: str):
    """
    Works even if the value is on the next line after the label.
    """
    pattern = re.compile(
        rf"{re.escape(tool_name)}\s+Predicted Stability Change\s*\(ΔΔG\)\s*:\s*(.*?)"
        rf"(?=(?:mCSM|SDM|DUET)\s+Predicted Stability Change\s*\(ΔΔG\)\s*:|Mutation:|$)",
        re.IGNORECASE | re.DOTALL
    )
    m = pattern.search(block)
    if not m:
        return None
    val = extract_first_float(m.group(1))
    return val

def extract_field(block: str, label: str):
    m = re.search(rf"{re.escape(label)}\s*:\s*([^\n]+)", block, re.IGNORECASE)
    return clean_text(m.group(1)) if m else ""

def normalize_triplet_aa(aa3: str) -> str:
    aa3 = aa3.strip().upper()
    return AA3_TO_1.get(aa3, "")

def build_record(block: str) -> dict:
    first_line = clean_text(block.splitlines()[0])
    hm = HEADER_RE.match(first_line)

    rec = {
        "uniprot_accession": "",
        "header_mutation": "",
        "mutation": "",
        "tested_mutation": "",
        "mutation_match": "",
        "original_aa": "",
        "position": pd.NA,
        "new_aa": "",
        "chain": "",
        "secondary_structure": "",
        "mcsm_ddg": pd.NA,
        "mcsm_interpretation": "",
        "sdm_ddg": pd.NA,
        "sdm_interpretation": "",
        "duet_ddg": pd.NA,
        "duet_interpretation": "",
        "raw_block": block,
    }

    if hm:
        rec["uniprot_accession"] = hm.group(1).strip()
        rec["header_mutation"] = hm.group(2).strip()

    rec["mcsm_ddg"] = extract_ddg(block, "mCSM")
    rec["sdm_ddg"] = extract_ddg(block, "SDM")
    rec["duet_ddg"] = extract_ddg(block, "DUET")

    rec["original_aa"] = extract_field(block, "Wild-type").upper()
    pos = extract_field(block, "Position")
    rec["position"] = int(pos) if pos.isdigit() else pd.NA
    rec["new_aa"] = extract_field(block, "Mutant-type").upper()
    rec["chain"] = extract_field(block, "Chain").upper()
    rec["secondary_structure"] = extract_field(block, "Secondary structure")

    old1 = normalize_triplet_aa(rec["original_aa"])
    new1 = normalize_triplet_aa(rec["new_aa"])

    if old1 and new1 and pd.notna(rec["position"]):
        rec["tested_mutation"] = f"{old1}{int(rec['position'])}{new1}"

    rec["mutation"] = rec["tested_mutation"] or rec["header_mutation"]

    if rec["header_mutation"] and rec["tested_mutation"]:
        rec["mutation_match"] = "match" if rec["header_mutation"] == rec["tested_mutation"] else "mismatch"

    for ddg_col, interp_col in [
        ("mcsm_ddg", "mcsm_interpretation"),
        ("sdm_ddg", "sdm_interpretation"),
        ("duet_ddg", "duet_interpretation"),
    ]:
        val = rec[ddg_col]
        if pd.notna(val):
            rec[interp_col] = "Stabilizing" if float(val) > 0 else "Destabilizing"

    return rec

def main():
    print(f"Loading: {DOCX_PATH}")
    text = load_docx_text(DOCX_PATH)
    blocks = split_into_blocks(text)

    records = [build_record(block) for block in blocks]
    df = pd.DataFrame(records)

    # remove clearly broken empty records
    df = df[(df["uniprot_accession"] != "") | (df["header_mutation"] != "")].copy()

    # remove exact duplicates by protein + tested mutation + duet value when possible
    dedup_cols = [c for c in ["uniprot_accession", "mutation", "duet_ddg"] if c in df.columns]
    if dedup_cols:
        df = df.drop_duplicates(subset=dedup_cols).copy()

    # sort cleanly
    if "position" in df.columns:
        df = df.sort_values(["uniprot_accession", "position"], na_position="last").reset_index(drop=True)

    ordered_cols = [
        "uniprot_accession",
        "header_mutation",
        "mutation",
        "tested_mutation",
        "mutation_match",
        "original_aa",
        "position",
        "new_aa",
        "chain",
        "secondary_structure",
        "mcsm_ddg",
        "mcsm_interpretation",
        "sdm_ddg",
        "sdm_interpretation",
        "duet_ddg",
        "duet_interpretation",
        "raw_block",
    ]
    df = df[[c for c in ordered_cols if c in df.columns]]

    print(f"Successfully parsed {len(df)} mutations")
    print("Columns extracted:", list(df.columns))
    print("\nFirst 5 rows:")
    print(df.head().to_string(index=False))

    print("\nFilled values:")
    for col in [
        "position", "new_aa", "chain", "secondary_structure",
        "mcsm_ddg", "sdm_ddg", "duet_ddg", "tested_mutation", "mutation_match"
    ]:
        if col in df.columns:
            print(f"  {col}: {df[col].notna().sum()}/{len(df)}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_PATH, index=False)
    print(f"\nSaved to: {OUT_PATH}")

if __name__ == "__main__":
    main()
