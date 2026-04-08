import docx
import pandas as pd
import re
from pathlib import Path

# ====================== LOAD DOCUMENT ======================
doc_path = "results/mutations/Duet_results.docx"
print(f"📂 Loading: {doc_path}")

doc = docx.Document(doc_path)

records = []
current = None
last_tool = None   # "mcsm", "sdm", or "duet"

for para in doc.paragraphs:
    text = para.text.strip()
    if not text:
        continue

    text_norm = text.replace('–', '-').replace('—', '-').replace('  ', ' ')

    # ==================== NEW MUTATION HEADER ====================
    if re.match(r'^[A-Z0-9a-z_]+', text_norm) and '-' in text_norm and len(text_norm) < 70:
        if current and 'uniprot_accession' in current:
            records.append(current.copy())
        
        current = {}
        parts = re.split(r'\s*[-–]\s*', text_norm, maxsplit=1)
        if len(parts) == 2:
            current['uniprot_accession'] = parts[0].strip()
            current['mutation'] = parts[1].strip()
        last_tool = None
        continue

    if current is None:
        continue

    # ==================== SET LAST TOOL (label lines) ====================
    if "mCSM Predicted Stability Change" in text:
        last_tool = "mcsm"
    elif "SDM Predicted Stability Change" in text:
        last_tool = "sdm"
    elif "DUET Predicted Stability Change" in text:
        last_tool = "duet"

    # ==================== EXTRACT VALUE ====================
    value_match = re.search(r'([+-]?\d*\.?\d+)', text)
    if value_match and last_tool:
        val = float(value_match.group(1))
        
        if last_tool == "mcsm":
            current['mcsm_ddg'] = val
        elif last_tool == "sdm":
            current['sdm_ddg'] = val
        elif last_tool == "duet":
            current['duet_ddg'] = val
        
        last_tool = None   # reset after using the value

    # ==================== MUTATION DETAILS ====================
    if "Wild-type" in text:
        m = re.search(r':\s*([A-Z]+)', text)
        if m:
            current['original_aa'] = m.group(1)

    if re.search(r'Position[:\s]', text):
        m = re.search(r'(\d+)', text)
        if m:
            current['position'] = int(m.group(1))

    if "Mutant-type" in text:
        m = re.search(r':\s*([A-Z]+)', text)
        if m:
            current['new_aa'] = m.group(1)

    if "Chain" in text:
        m = re.search(r':\s*([A-Z])', text)
        if m:
            current['chain'] = m.group(1)

    if "Secondary structure" in text:
        current['secondary_structure'] = text.split(":", 1)[1].strip() if ":" in text else text.split()[-1]

# Save last record
if current and 'uniprot_accession' in current:
    records.append(current.copy())

# ====================== DATAFRAME ======================
df = pd.DataFrame(records)

# Add interpretation columns
for col in ['mcsm_ddg', 'sdm_ddg', 'duet_ddg']:
    if col in df.columns:
        interp = col.replace('_ddg', '_interpretation')
        df[interp] = df[col].apply(lambda x: 'Stabilizing' if x > 0 else 'Destabilizing')

# Clean up
df = df.sort_values(by=['uniprot_accession', 'position']).reset_index(drop=True)

print(f"✅ Successfully parsed {len(df)} mutations")
print("Columns extracted:", list(df.columns))
print("\nFirst 5 rows:")
print(df.head())

# Quick summary of filled columns
print("\nFilled values:")
for col in ['mcsm_ddg', 'sdm_ddg', 'duet_ddg']:
    if col in df.columns:
        filled = df[col].notna().sum()
        print(f"  {col}: {filled}/{len(df)}")

# ====================== SAVE ======================
output_path = Path("results/mutations/duet_results_manual.csv")
output_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_path, index=False)

print(f"\n📁 Saved to: {output_path}")
