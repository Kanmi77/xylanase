#!/usr/bin/env bash
set -euo pipefail

BASE=~/xylanase/xylanase
INDIR="$BASE/results/alignments"
OUTDIR="$BASE/results/trees"
LOGDIR="$BASE/results/logs/fasttree"

mkdir -p "$OUTDIR" "$LOGDIR"

for aln in "$INDIR"/*_aligned.fasta
do
  [ -s "$aln" ] || continue
  base=$(basename "$aln" _aligned.fasta)
  FastTree -lg "$aln" > "$OUTDIR/${base}.nwk" 2> "$LOGDIR/${base}_fasttree.log"
  echo "[DONE] $base"
done
