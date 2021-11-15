#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

if [ ! -f JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt ]; then
  wget "http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt"
fi

if [ ! -d JASPAR2020_CORE_plants_non-redundant_pfms_meme ]; then
  wget "http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_plants_non-redundant_pfms_meme.zip"
  unzip -d JASPAR2020_CORE_plants_non-redundant_pfms_meme JASPAR2020_CORE_plants_non-redundant_pfms_meme.zip
  rm JASPAR2020_CORE_plants_non-redundant_pfms_meme.zip
fi

if [ ! -f clusters.tsv ]; then
  wget "http://jaspar.genereg.net/static/clustering/JASPAR_2020_clusters/plants/interactive_trees/JASPAR_2020_matrix_clustering_plants_archive.zip"
  unzip -p JASPAR_2020_matrix_clustering_plants_archive.zip interactive_trees/JASPAR_2020_matrix_clustering_plants_tables/clusters.tab | sed -E 's/JASPAR_2020_plants_m[0-9]+_(MA[0-9]+)_([0-9]+)/\1.\2/g' > clusters.tsv
  rm JASPAR_2020_matrix_clustering_plants_archive.zip
fi

chmod -wx clusters.tsv JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt JASPAR2020_CORE_plants_non-redundant_pfms_meme/*

