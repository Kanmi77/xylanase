#!/usr/bin/env bash
set -euo pipefail

BASE=~/xylanase/xylanase
INDIR="$BASE/results/docking/cleaned_receptors"
OUTDIR="$BASE/results/docking/receptors_pdbqt"
LOGDIR="$BASE/results/docking/logs"

mkdir -p "$OUTDIR" "$LOGDIR"

for pdb in "$INDIR"/*.pdb
do
  [ -s "$pdb" ] || continue
  base=$(basename "$pdb" .pdb)

  mk_prepare_receptor.py \
    --read_pdb "$pdb" \
    -o "$OUTDIR/${base}" \
    -p \
    > "$LOGDIR/${base}_prepare_receptor.log" 2>&1

  echo "[DONE] $base"
done
