#!/usr/bin/env bash
set -euo pipefail

IDX=""
LIST=""
OUT=""
THREADS=4

usage() {
  echo "Usage: $0 -x <bt2_index_prefix> -l <fastq_list.txt> -o <outdir> [-t <threads>]"
  echo "  fastq_list.txt: two columns per line: R1 <tab or space> R2"
  exit 1
}

while getopts "x:l:o:t:" opt; do
  case $opt in
    x) IDX="$OPTARG" ;;
    l) LIST="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -n "$IDX" && -n "$LIST" && -n "$OUT" ]] || usage
[[ -r "$LIST" ]] || { echo "ERROR: Cannot read list file: $LIST"; exit 1; }

if [[ ! -f "${IDX}.1.bt2" && ! -f "${IDX}.1.bt2l" ]]; then
  echo "ERROR: Bowtie2 index not found for prefix: $IDX"
  exit 1
fi

mkdir -p "$OUT"

echo "Index prefix : $IDX"
echo "List file    : $LIST"
echo "Output dir   : $OUT"
echo "Threads      : $THREADS"
echo

while IFS=$' \t' read -r R1 R2 REST || [[ -n "${R1:-}" ]]; do
  R1="${R1//$'\r'/}"
  R2="${R2//$'\r'/}"
  [[ -z "${R1:-}" ]] && continue
  [[ "${R1:0:1}" == "#" ]] && continue
  [[ -z "${R2:-}" ]] && { echo "WARNING: missing R2 for $R1"; continue; }

  [[ -f "$R1" ]] || { echo "WARNING: R1 not found: $R1"; continue; }
  [[ -f "$R2" ]] || { echo "WARNING: R2 not found: $R2"; continue; }

  base=$(basename "$R1")
  SAMPLE="${base%.fastq.gz}"
  SAMPLE="${SAMPLE%.fq.gz}"
  SAMPLE="${SAMPLE%_R1_001}"
  SAMPLE="${SAMPLE%_R1}"
  SAMPLE="${SAMPLE%_1}"

  echo "▶ Processing: $SAMPLE"

  bowtie2 -x "$IDX" \
    -1 "$R1" -2 "$R2" \
    -p "$THREADS" \
    --very-sensitive \
    --un-conc-gz "$OUT/${SAMPLE}.hostRemoved.%.fastq.gz" \
    -S /dev/null \
    2> "$OUT/${SAMPLE}.bowtie2.log"

  mv "$OUT/${SAMPLE}.hostRemoved.1.fastq.gz" "$OUT/${SAMPLE}_R1.clean.fastq.gz"
  mv "$OUT/${SAMPLE}.hostRemoved.2.fastq.gz" "$OUT/${SAMPLE}_R2.clean.fastq.gz"

  echo "  ✓ Wrote: $OUT/${SAMPLE}_R1.clean.fastq.gz"
  echo "  ✓ Wrote: $OUT/${SAMPLE}_R2.clean.fastq.gz"
  echo
done < "$LIST"

echo "Done."
