#!/usr/bin/env bash
set -euo pipefail

OUTDIR=~/xylanase/data/refseq/raw
mkdir -p "$OUTDIR"

echo "Prepare RefSeq acquisition here." > "$OUTDIR/README_refseq.txt"
echo "This folder will store RefSeq protein accessions, fasta files, and metadata." >> "$OUTDIR/README_refseq.txt"

echo "Created RefSeq placeholder at $OUTDIR/README_refseq.txt"

