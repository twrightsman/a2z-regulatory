#!/bin/bash

set -euo pipefail

if [ ! -d "tmp" ]; then
  mkdir tmp
fi

WINDOW_SIZE="${1:-}"
INTERVAL_FILE="${2:-}"
IDX="${3:-}"

if [ -z "$WINDOW_SIZE" ] || [ -z "$INTERVAL_FILE" ] || [ -z "$IDX" ]; then
  echo "Usage: $0 WINDOW_SIZE peaks.bed genome.fa.fai" >&2
  exit 1
fi

INTERVAL_FILE_BASE=$(basename $INTERVAL_FILE)
SPECIES="${INTERVAL_FILE_BASE%%.*}"

TMP=$(mktemp -d)
echo "Using temporary directory $TMP" >&2

if (( (WINDOW_SIZE % 2) != 0 )) || (( WINDOW_SIZE <= 1 )); then
  echo "Window size needs to be even and 2bp or more" >&2
  return 1
else
  echo "Using window size of ${WINDOW_SIZE}bp" >&2
fi

# grab seqid, start, and end of positive peaks
# expand peaks to WINDOW_SIZE, from midpoint
awk 'NR==FNR{A[$1];next}($1 in A)' "$IDX" "${INTERVAL_FILE}" | \
  cut -f1-3 > "${TMP}/${SPECIES}.pos.filt.bed"
awk -v OFS=$'\t' '{mid = int(($2+$3)/2); print $1, mid, mid + 1}' < "${TMP}/${SPECIES}.pos.filt.bed" | \
  bedtools slop -l $(( WINDOW_SIZE / 2 )) -r $(( (WINDOW_SIZE / 2) - 1 )) -g "$IDX" -i - | \
  awk '($3 - $2) == '$WINDOW_SIZE > "${TMP}/${SPECIES}.pos.${WINDOW_SIZE}.bed"
bedtools sort -i "${TMP}/${SPECIES}.pos.${WINDOW_SIZE}.bed" -g "$IDX" > "${TMP}/${SPECIES}.pos.bed"

# get the complement of the resized shorter and unresized longer open regions and make non-overlapping windows of WINDOW_SIZE
# filter out sequences that don't have at least one open region
awk '($3 - $2) > '$WINDOW_SIZE < "${TMP}/${SPECIES}.pos.filt.bed" > "${TMP}/${SPECIES}.pos.large.bed"
bedtools sort -i <(cat "${TMP}/${SPECIES}.pos.${WINDOW_SIZE}.bed" "${TMP}/${SPECIES}.pos.large.bed") -g "$IDX" | \
  bedtools merge -i - | \
  bedtools complement -i - -g "$IDX" | \
  bedtools makewindows -b - -w $WINDOW_SIZE | \
  awk '($3 - $2) == '$WINDOW_SIZE | \
  cut -f1-3 | \
  awk 'NR==FNR{A[$1];next}($1 in A)' <(cut -f 1 "${INTERVAL_FILE}" | sort | uniq) - > "${TMP}/${SPECIES}.neg.bed"

# add targets/labels, trim windows on edges of chromosomes (for safety during random shifting), and combine into one file
cat <(sed 's/$/\t1/' "${TMP}/${SPECIES}.pos.bed") <(sed 's/$/\t0/' "${TMP}/${SPECIES}.neg.bed") | \
  bedtools intersect -wa -f 1 -a - -b <(awk -v OFS=$'\t' '{if ($2 >= 1000) {print $1, 100, $2 - 100}}' "$IDX") | \
  bedtools sort -i - > "tmp/${SPECIES}.unlabeled.bed"

