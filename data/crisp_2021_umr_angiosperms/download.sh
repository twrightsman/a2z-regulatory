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

# References and annotations
mkdir -p {genomes,annotations,intervals}

## Sorghum bicolor (v3.0)
jgi_download "Sbicolor" "Sbicolor_454_v3.1.1.gene.gff3.gz" "annotations/Sorghum_bicolor.gff3"

## Zea mays (AGPv4)
if [ ! -f annotations/Zea_mays.gff3 ]; then
  curl "https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz" | gzip -cd | sed 's/^Chr//' | grep -v '^Mt' > annotations/Zea_mays.gff3
fi

## Vitis vinifera (Genoscope.12X, Phytozome 10)
if [ ! -f genomes/Vitis_vinifera.fa ]; then
  wget "https://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/assembly/goldenpath/unmasked/chr.fa" -O genomes/Vitis_vinifera.fa
fi
if [ ! -f annotations/Vitis_vinifera.gff3 ]; then
  wget "https://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/annotation/Vitis_vinifera_annotation.gff.gz" -O - | gzip -cd > annotations/Vitis_vinifera.gff3
  samtools faidx genomes/Vitis_vinifera.fa
fi

# Make data read-only for safety
chmod -wx {annotations,genomes}/*

# Link UMRs
ln --force --symbolic ../raw/Arabidopsis_thaliana_II/Athaliana_cov_3_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed intervals/Arabidopsis_thaliana.bed
ln --force --symbolic ../raw/Eutrema_salsugineum_II/Eutrema_salsugineum_cov_3_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed intervals/Eutrema_salsugineum.bed
ln --force --symbolic ../raw/Glycine_max_II/Gmax_LD_cov_3_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed intervals/Glycine_max.bed
ln --force --symbolic ../raw/Populus_trichocarpa_II/Populus_trichocarpa_cov_3_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed intervals/Populus_trichocarpa.bed
ln --force --symbolic ../raw/Vitis_vinifera_II/Vvinifera_cov_3_sites_2_MR_0.4_UMR_0.1_UMRs_6col.bed intervals/Vitis_vinifera.bed

ln --force --symbolic ../../crisp_2020_umr_grasses/intervals/Brachypodium_distachyon.bed intervals/
ln --force --symbolic ../../crisp_2020_umr_grasses/intervals/Hordeum_vulgare.bed intervals/
ln --force --symbolic ../../crisp_2020_umr_grasses/intervals/Oryza_sativa.bed intervals/
ln --force --symbolic ../../crisp_2020_umr_grasses/intervals/Sorghum_bicolor.bed intervals/
ln --force --symbolic ../../crisp_2020_umr_grasses/intervals/Zea_mays.bed intervals/

# Link annotations
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Arabidopsis_thaliana.gff3 annotations/
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Brachypodium_distachyon.gff3 annotations/
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Eutrema_salsugineum.gff3 annotations/
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Glycine_max.gff3 annotations/
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Hordeum_vulgare.gff3 annotations/
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Oryza_sativa.gff3 annotations/
ln --force --symbolic ../../lu_2019_ATACseq/annotations/Populus_trichocarpa.gff3 annotations/

# Link genomes
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Arabidopsis_thaliana.fa genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Arabidopsis_thaliana.fa.fai genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Brachypodium_distachyon.fa genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Brachypodium_distachyon.fa.fai genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Eutrema_salsugineum.fa genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Eutrema_salsugineum.fa.fai genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Glycine_max.fa genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Glycine_max.fa.fai genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Hordeum_vulgare.fa genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Hordeum_vulgare.fa.fai genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Oryza_sativa.fa genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Oryza_sativa.fa.fai genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Populus_trichocarpa.fa genomes/
ln --force --symbolic ../../lu_2019_ATACseq/genomes/Populus_trichocarpa.fa.fai genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Sorghum_bicolor.fa genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Sorghum_bicolor.fa.fai genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Zea_mays.fa genomes/
ln --force --symbolic ../../crisp_2020_umr_grasses/genomes/Zea_mays.fa.fai genomes/

# Delete session credentials once finished
rm --force cookies.txt
