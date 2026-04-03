#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
from Bio import AlignIO
import pandas as pd

BASE = Path.home() / "xylanase" / "xylanase"
ALN = BASE / "results" / "mutations" / "thermostable_labeled_subset_aligned.fasta"
OUT = BASE / "results" / "mutations" / "thermostable_consensus_positions.csv"

alignment = AlignIO.read(ALN, "fasta")
records = []

for i in range(alignment.get_alignment_length()):
    col = alignment[:, i]
    aa = [x for x in col if x != "-"]
    if not aa:
        continue
    counts = Counter(aa)
    consensus_aa, freq = counts.most_common(1)[0]
    records.append({
        "alignment_position": i + 1,
        "consensus_aa": consensus_aa,
        "consensus_count": freq,
        "non_gap_count": len(aa),
        "consensus_fraction": freq / len(aa),
        "column_residues": "".join(col),
    })

pd.DataFrame(records).to_csv(OUT, index=False)
print(f"Saved: {OUT}")
