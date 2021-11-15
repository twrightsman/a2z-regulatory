#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

if [ ! -f lin_all.peaks.final.bed ] || [ ! -f GM.peaks.final.bed ]; then
  wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE66386&format=file&file=GSE66386_all_nucleosome_calling.tar.gz' -O - | gzip -cd | tar xvf -
  rm osmotic.peaks.final.bed pombe.peaks.final.bed
fi

if [ ! -f yeast.renamed.fa ]; then
  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_report.txt' -O - | grep -v '^#' | cut -f 1,7 > yeast_name2refseq.tsv
  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz' -O yeast.ncbi.fa.gz
  gzip -cd yeast.ncbi.fa.gz | awk '{if (NR==FNR){A[$2] = $1;next} else {if ($0 ~ /^>/){print ">" A[substr($1, 2)]} else {print $0}}}' yeast_name2refseq.tsv - | sed -E 's/^>([IVX]+)/>chr\1/' > yeast.renamed.fa
  rm yeast_name2refseq.tsv yeast.ncbi.fa.gz
  samtools faidx yeast.renamed.fa
fi

if [ ! -f hg19.fa ]; then
  for i in {1..22}; do
    wget 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chr'$i'.fa.gz'
  done
  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chrX.fa.gz'
  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chrY.fa.gz'
  zcat chr*.fa.gz | sed -E 's/^>.*chromosome ([0-9XY]+).*/>chr\1/' > hg19.fa
  rm chr*.fa.gz
  samtools faidx hg19.fa
fi

ln --symbolic --force lin_all.peaks.final.bed yeast.open.bed

chmod -wx *.bed *.fa*
