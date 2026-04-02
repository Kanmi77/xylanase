#!/usr/bin/env bash
set -euo pipefail

OUTDIR=~/xylanase/data/cazy/raw
mkdir -p "$OUTDIR"

echo "Prepare CAZy acquisition here." > "$OUTDIR/README_cazy.txt"
echo "This folder will store GH10/GH11 xylanase records and extracted metadata." >> "$OUTDIR/README_cazy.txt"

echo "Created CAZy placeholder at $OUTDIR/README_cazy.txt"
