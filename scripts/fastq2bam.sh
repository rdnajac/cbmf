#!/bin/bash
# run this script in the folder to process
set -euo pipefail
set -x
THREADS=$(nproc)
SUFFIX="_001.fastq.gz"
# REFERENCE_GENOME="$HUMANREF"
REFERENCE_GENOME="$MOUSEREF"

function get_fastq_pairs {
  local dir=${1:-$(pwd)}
  # trim the trailing slash
  if [[ "${dir}" == */ ]]; then
    dir="${dir%/}"
  fi
  declare -a files=()
  for r1 in "${dir}"/*_R1_001.fastq.gz; do
    local r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    [[ -f "$r2" ]] && files+=("${r1%_R1_001.fastq.gz}") || echo "error: ${r1##*/}"
  done
  printf "%s\n" "${files[@]}"
}

fastq2bam() {
  local ref=$1
  local id=$2

bowtie2 -t -p "${THREADS}" --mm -x "${ref}" -1 "${id}_R1_001.fastq.gz" -2 "${id}_R2_001.fastq.gz" 2>> "${id}.log" | samtools sort -@ "${THREADS}" -o "${id}.bam" 2>> "${id}.log"
}

declare -a RA=$(get_fastq_pairs ra/)

main() {
  for sample in ${RA[@]}; do
    fastq2bam "${REFERENCE_GENOME}" "${sample}"
  done
}

main
