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

# References
mkdir -p {genomes,intervals}

## Brachypodium distachyon (v3.1)
ln --symbolic --force ../../lu_2019_ATACseq/genomes/Brachypodium_distachyon.fa genomes/Brachypodium_distachyon.fa
ln --symbolic --force ../../lu_2019_ATACseq/genomes/Brachypodium_distachyon.fa.fai genomes/Brachypodium_distachyon.fa.fai

## Oryza sativa (Rice, v7.0)
ln --symbolic --force ../../lu_2019_ATACseq/genomes/Oryza_sativa.fa genomes/Oryza_sativa.fa
ln --symbolic --force ../../lu_2019_ATACseq/genomes/Oryza_sativa.fa.fai genomes/Oryza_sativa.fa.fai

## Sorghum bicolor (v3.0)
jgi_download "Sbicolor" "Sbicolor_454_v3.0.1.fa.gz" "genomes/Sorghum_bicolor.fa"
if [ ! -f genomes/Sorghum_bicolor.fa.fai ]; then
  samtools faidx genomes/Sorghum_bicolor.fa
fi

## Zea mays (Maize/Corn, AGPv4)
if [ ! -f genomes/Zea_mays.fa ]; then
  curl "https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz" | gzip -cd | sed 's/^>Chr/>/' > genomes/Zea_mays.fa
  samtools faidx genomes/Zea_mays.fa
fi

## Hordeum vulgare L. (Barley, r42)
ln --symbolic --force ../../lu_2019_ATACseq/genomes/Hordeum_vulgare.fa genomes/Hordeum_vulgare.fa
ln --symbolic --force ../../lu_2019_ATACseq/genomes/Hordeum_vulgare.fa.fai genomes/Hordeum_vulgare.fa.fai

# UMRs
if [ ! -f intervals/Zea_mays.bed ]; then
  wget https://www.pnas.org/highwire/filestream/947726/field_highwire_adjunct_files/3/pnas.2010250117.sd03.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > intervals/Zea_mays.bed
fi
if [ ! -f intervals/Brachypodium_distachyon.bed ]; then
  wget https://www.pnas.org/highwire/filestream/947726/field_highwire_adjunct_files/10/pnas.2010250117.sd10.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > intervals/Brachypodium_distachyon.bed
fi
if [ ! -f intervals/Hordeum_vulgare.bed ]; then
  wget https://www.pnas.org/highwire/filestream/947726/field_highwire_adjunct_files/7/pnas.2010250117.sd07.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > intervals/Hordeum_vulgare.bed
fi
if [ ! -f intervals/Oryza_sativa.bed ]; then
  wget https://www.pnas.org/highwire/filestream/947726/field_highwire_adjunct_files/9/pnas.2010250117.sd09.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > intervals/Oryza_sativa.bed
fi
if [ ! -f intervals/Sorghum_bicolor.bed ]; then
  wget https://www.pnas.org/highwire/filestream/947726/field_highwire_adjunct_files/8/pnas.2010250117.sd08.csv -O - | tail -n+2 | tr ',' $'\t' | cut -f1-3 > intervals/Sorghum_bicolor.bed
fi

# Make data read-only for safety
chmod -wx {intervals,genomes}/*

# Delete session credentials once finished
rm --force cookies.txt
