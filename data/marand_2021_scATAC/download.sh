#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

if [ ! -f GSE155178_maize_scATAC_atlas_ACR_celltype_CPM.txt ]; then
  wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155178&format=file&file=GSE155178%5Fmaize%5FscATAC%5Fatlas%5FACR%5Fcelltype%5FCPM%2Etxt%2Egz" -O GSE155178_maize_scATAC_atlas_ACR_celltype_CPM.txt.gz
  gzip -d GSE155178_maize_scATAC_atlas_ACR_celltype_CPM.txt.gz
fi

if [ ! -f AGPv4.fa ]; then
  wget "ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz" -O - | \
    gzip -cd | \
    sed -r 's/^>([0-9]+)/>chr\1/' > AGPv4.fa
  samtools faidx AGPv4.fa
fi

chmod -wx GSE155178_maize_scATAC_atlas_ACR_celltype_CPM.txt AGPv4.fa*
