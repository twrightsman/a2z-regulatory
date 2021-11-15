#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

if [ ! -f v4.repeats.renamed.bed ]; then
  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/4577/102/GCF_000005005.2_B73_RefGen_v4/GCF_000005005.2_B73_RefGen_v4_rm.out.gz' -O - | \
    gzip -cd | \
    tail -n+4 | \
    tr -s ' ' | \
    sed 's/^ //' | \
    cut -d ' ' -f 5-7 --output-delimiter=$'\t' > v4.repeats.bed
  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/4577/102/GCF_000005005.2_B73_RefGen_v4/GCF_000005005.2_B73_RefGen_v4_assembly_report.txt' -O - | \
    grep -v '^#' | \
    cut -f 1,7 > name2refseq.tsv
  awk -v 'OFS=\t' '{ if (NR==FNR){A[$2] = $1;next} else {print A[$1], $2, $3} }' name2refseq.tsv v4.repeats.bed | sed 's/^chr//' > v4.repeats.renamed.bed
  rm v4.repeats.bed name2refseq.tsv
fi

chmod -wx v4.repeats.renamed.bed
