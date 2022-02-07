#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# log in to JGI
COOKIES="cookies.txt"
function jgi_login() {
  if [ ! -f "$COOKIES" ]; then
    read -p "JGI Phytozome Username: " JGI_USER
    read -s -p "JGI Phytozome Password: " JGI_PASSWORD
    echo
    curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode "login=${JGI_USER}" --data-urlencode "password=${JGI_PASSWORD}" --cookie-jar "$COOKIES" > /dev/null 2>&1
  fi
}

function jgi_download() {
  local portal="$1"
  local filename="$2"
  local out_path="$3"
  downloaded=0

  if [ ! -f "$out_path" ]; then
    jgi_login
    echo "Downloading $filename from JGI portal $portal"
    url=$(curl -b "$COOKIES" "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=${portal}" 2> /dev/null | xmllint --xpath "//file[@filename='${filename}']/@url" /dev/stdin | sed -E 's/ *url="(.*)"$/\1/' | sed 's/&amp;/\&/')
    curl --cookie "$COOKIES" "https://genome.jgi.doe.gov${url}" | gzip --stdout --decompress > "$out_path"
    downloaded=1
  else
    echo "$out_path already exists, skipping"
  fi
}

# ACRs
PEAKS="peaks"
if [ ! -d "$PEAKS" ]; then
  mkdir "$PEAKS"
  wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE128434&format=file&file=GSE128434_BED_files.tar.gz' --output-document=GSE128434_BED_files.tar.gz
  tar --extract --file=GSE128434_BED_files.tar.gz
  rm GSE128434_BED_files.tar.gz
  gzip --stdout --decompress Arabidopsis_7days_leaf_ACRs.bed.gz > "${PEAKS}/Arabidopsis_thaliana.bed"
  gzip --stdout --decompress Asparagus_10days_leaf_ACRs.bed.gz > "${PEAKS}/Asparagus_officinalis.bed"
  gzip --stdout --decompress Barley_7days_leaf_ACRs.bed.gz > "${PEAKS}/Hordeum_vulgare.bed"
  gzip --stdout --decompress Brachypodium_7days_leaf_ACRs.bed.gz > "${PEAKS}/Brachypodium_distachyon.bed"
  gzip --stdout --decompress Eutrema_10days_leaf_ACRs.bed.gz > "${PEAKS}/Eutrema_salsugineum.bed"
  gzip --stdout --decompress Maize_7days_leaf_ACRs.bed.gz > "${PEAKS}/Zea_mays.bed"
  gzip --stdout --decompress Phaseolus_10days_leaf_ACRs.bed.gz > "${PEAKS}/Phaseolus_vulgaris.bed"
  gzip --stdout --decompress Populus_leaf_ACRs.bed.gz > "${PEAKS}/Populus_trichocarpa.bed"
  gzip --stdout --decompress Rice_7days_leaf_ACRs.bed.gz > "${PEAKS}/Oryza_sativa.bed"
  gzip --stdout --decompress Setaria_7days_leaf_ACRs.bed.gz > "${PEAKS}/Setaria_viridis.bed"
  gzip --stdout --decompress Sorghum_7days_leaf_ACRs.bed.gz > "${PEAKS}/Sorghum_bicolor.bed"
  gzip --stdout --decompress Soybean_10days_leaf_ACRs.bed.gz > "${PEAKS}/Glycine_max.bed"
  rm *.bed.gz
fi


# References and annotations
mkdir -p {genomes,peaks,annotations}

function download_decmp() {
  local url="$1"
  local out_path="$2"
  downloaded=0

  if [ ! -f "$out_path" ]; then
    echo "Downloading $url"
    curl "$url" | gzip --stdout --decompress > "$out_path"
    downloaded=1
  else
    echo "$out_path already exists, skipping"
  fi
}

## Arabidopsis thaliana (TAIR10)
jgi_download "Athaliana" "Athaliana_447_TAIR10.fa.gz" "genomes/Arabidopsis_thaliana.fa"
if [ $downloaded -eq 1 ]; then
  # rename chromosomes to lowercase to match BED
  sed --in-place 's/^>Chr/>chr/' "genomes/Arabidopsis_thaliana.fa"
fi
jgi_download "Athaliana" "Athaliana_167_TAIR10.gene.gff3.gz" "annotations/Arabidopsis_thaliana.gff3"
if [ $downloaded -eq 1 ]; then
  sed --in-place 's/^Chr/chr/' "annotations/Arabidopsis_thaliana.gff3"
fi

## Eutrema salsugineum (v1.0)
jgi_download "Esalsugineum" "Esalsugineum_173_v1.fa.gz" "genomes/Eutrema_salsugineum.fa"
jgi_download "Esalsugineum" "Esalsugineum_173_v1.0.gene.gff3.gz" "annotations/Eutrema_salsugineum.gff3"

## Phaseolus vulgaris (Common bean, v1.0)
jgi_download "Pvulgaris" "Pvulgaris_218_v1.0.fa.gz" "genomes/Phaseolus_vulgaris.fa"
jgi_download "Pvulgaris" "Pvulgaris_218_v1.0.gene.gff3.gz" "annotations/Phaseolus_vulgaris.gff3"

## Glycine max (Soybean, Wm82.a2.v1)
jgi_download "Gmax" "Gmax_275_v2.0.fa.gz" "genomes/Glycine_max.fa"
if [ $downloaded -eq 1 ]; then
  sed --in-place --regexp-extended 's/^>Chr0?([0-9]+)/>chr\1/' "genomes/Glycine_max.fa"
fi
jgi_download "Gmax" "Gmax_275_Wm82.a2.v1.gene.gff3.gz" "annotations/Glycine_max.gff3"
if [ $downloaded -eq 1 ]; then
  sed --in-place --regexp-extended 's/^Chr0?([0-9]+)/chr\1/' "annotations/Glycine_max.gff3"
fi

## Brachypodium distachyon (v3.0)
jgi_download "Bdistachyon" "Bdistachyon_307_v3.0.fa.gz" "genomes/Brachypodium_distachyon.fa"
jgi_download "Bdistachyon" "Bdistachyon_307_v3.0.gene.gff3.gz" "annotations/Brachypodium_distachyon.gff3"

## Oryza sativa (Rice, v7.0)
jgi_download "Osativa" "Osativa_323_v7.0.fa.gz" "genomes/Oryza_sativa.fa"
jgi_download "Osativa" "Osativa_323_v7.0.gene.gff3.gz" "annotations/Oryza_sativa.gff3"

## Setaria viridis (v1.0)
jgi_download "Sviridis" "Sviridis_311_v1.0.fa.gz" "genomes/Setaria_viridis.fa"
jgi_download "Sviridis" "Sviridis_311_v1.1.gene.gff3.gz" "annotations/Setaria_viridis.gff3"

## Populus trichocarpa (Poplar, v3.0)
jgi_download "Ptrichocarpa" "Ptrichocarpa_210_v3.0.fa.gz" "genomes/Populus_trichocarpa.fa"
jgi_download "Ptrichocarpa" "Ptrichocarpa_210_v3.0.gene.gff3.gz" "annotations/Populus_trichocarpa.gff3"

## Sorghum bicolor (v3.0)
jgi_download "Sbicolor" "Sbicolor_313_v3.0.fa.gz" "genomes/Sorghum_bicolor.fa"
jgi_download "Sbicolor" "Sbicolor_313_v3.1.gene.gff3.gz" "annotations/Sorghum_bicolor.gff3"

## Zea mays (Maize/Corn, AGPv4.38)
download_decmp "ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz" "genomes/Zea_mays.fa"
download_decmp "ftp://ftp.ensemblgenomes.org/pub/plants/release-38/gff3/zea_mays/Zea_mays.AGPv4.38.gff3.gz" "annotations/Zea_mays.gff3"

## Hordeum vulgare L. (Barley, IBSC_v2)
download_decmp "ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa.gz" "genomes/Hordeum_vulgare.fa"
download_decmp "ftp://ftp.ensemblgenomes.org/pub/plants/release-42/gff3/hordeum_vulgare/Hordeum_vulgare.IBSC_v2.42.gff3.gz" "annotations/Hordeum_vulgare.gff3"

## Asparagus officinalis (Asparagus, v1.1)
download_decmp "http://asparagus.uga.edu/genome_files/AsparagusCHR_V1.1.fsa.gz" "genomes/Asparagus_officinalis.fa"
download_decmp "http://asparagus.uga.edu/genome_files/AsparagusCHR_V1.1.Final.gff3.gz" "annotations/Asparagus_officinalis.gff3"

for genome in genomes/*.fa; do
  index="${genome}.fai"
  if [ ! -f "$index" ] || [ "$genome" -nt "$index" ]; then
    echo "Indexing $genome"
    rm --force "$index"
    samtools faidx $genome
  fi
done

# Make data read-only for safety
chmod -wx {annotations,genomes,peaks}/* ACR_column_header.txt

# Delete session credentials once finished
rm --force cookies.txt
