#!/usr/bin/env bash
set -euo pipefail

BASE=~/xylanase/xylanase
INDIR="$BASE/data/phylogeny"
OUTDIR="$BASE/results/alignments"
LOGDIR="$BASE/results/logs/mafft"

mkdir -p "$OUTDIR" "$LOGDIR"

for fasta in \
  "$INDIR/all_filtered_xylanases.fasta" \
  "$INDIR/bacterial_GH10.fasta" \
  "$INDIR/bacterial_GH11.fasta" \
  "$INDIR/fungal_GH10.fasta" \
  "$INDIR/fungal_GH11.fasta"
do
  [ -s "$fasta" ] || continue
  base=$(basename "$fasta" .fasta)
  mafft --auto "$fasta" > "$OUTDIR/${base}_aligned.fasta" 2> "$LOGDIR/${base}_mafft.log"
  echo "[DONE] $base"
done
