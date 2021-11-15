#!/bin/bash

set -euo pipefail

IN_DIR="$1"
MEME_FILE="$2"

tail -n+2 "${IN_DIR}/high_effect_kmers.tsv" | \
  awk '{print ">kmer" NR "\n" $1}' > "${IN_DIR}/high_effect_kmers.fa"

fimo \
  --max-strand \
  --oc "${IN_DIR}/high_effect_fimo" \
  --qv-thresh \
  --thresh "0.05" \
  ../../data/JASPAR/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt \
  "${IN_DIR}/high_effect_kmers.fa"
