#!/usr/bin/env bash
set -euo pipefail

mkdir -p ~/xylanase/{data/{uniprot,refseq,cazy,raw,interim,curated,metadata},scripts/{00_setup,01_data_acquisition,02_curation,03_annotation,04_homology_models},models/{swiss_model,templates,qa},results/{logs,reports},docs}

mkdir -p ~/xylanase/data/uniprot/{raw,processed}
mkdir -p ~/xylanase/data/refseq/{raw,processed}
mkdir -p ~/xylanase/data/cazy/{raw,processed}
mkdir -p ~/xylanase/models/swiss_model/{inputs,outputs}
mkdir -p ~/xylanase/models/templates/{pdb,metadata}
mkdir -p ~/xylanase/results/logs/{downloads,curation,annotation,modeling}

echo "Project structure created at ~/xylanase"
find ~/xylanase -maxdepth 4 -type d | sort
