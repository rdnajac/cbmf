#!/bin/bash
#================================================================
## Description: Fastq to aligned and BAM or CRAM alignment pipeline
## Author: Ryan D. Najac
## Last modified: 2024-05-03
#================================================================
set -euo pipefail   # enforce strict error handling
set -x

#======================= GLOBALS ================================
# define these after you download the reference genomes
readonly HUMANREF="/home/ubuntu/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"
readonly MOUSEREF="/home/ubuntu/genomes/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index"

# use the maximum number of threads available
THREADS=$(nproc)

#======================= UTILITY =================================

# Helper functions to print colored messages
okay() { printf "\033[92m%s\033[0m\n" "${1}" >&2; }
info() { printf "\033[94m%s\033[0m\n" "${1}" >&2; }
warn() { printf "\033[91m%s\033[0m\n" "${1}" >&2; }

# Quit with an error message
bail() { warn "${1}"; exit "${2:-1}"; }

# Print usage and exit
usage() { printf "%s\n" "$USAGE"; exit 1; }

#======================= MAIN FUNCTIONALITY =====================
get_fastq_pairs() {
  local dir=${1:-$(pwd)}
  local suffix="_001.fastq.gz"
  declare -a files=()

  while IFS= read -r r1; do
    sample_id=$(basename "$r1" "_R1$suffix")
    local r2="$dir/${sample_id}_R2$suffix"
    [[ -f "$r2" ]] && files+=("${dir}/${sample_id}")
  done < <(find "$dir" -maxdepth 1 -type f -regex ".*_R[12]$suffix")

  printf "%s\n" "${files[@]}"
}

fastq2bam() {
  local id="$1"
  local r1="${id}_R1_001.fastq.gz"
  local r2="${id}_R2_001.fastq.gz"
  local log_file="${id}.log"

  log() { echo "$(date +%H%M%S) $*" >> "$log_file"; }

  log "start $id"
  {
    # Align the reads, sort by name, fix mate information, sort again, and mark duplicates
    bowtie2 --time --threads "$THREADS" --mm -x "$MOUSEREF" -1 "$r1" -2 "$r2" \
      | samtools sort -n -@ "$THREADS" - \
      | samtools fixmate -m -@ "$THREADS" - - \
      | samtools sort -@ "$THREADS" - \
      | samtools markdup -@ "$THREADS" - "${id}.bam"
  } 2> "$log_file"
  log "finish $id"
}

#======================= MAIN ====================================
#[[ s -d "$1" ]] && bail "Error: Directory '$1' not found."

declare -a fastq_files=()

# kind of hacky, but this will work for now
fastq_files=($(get_fastq_pairs $1))

for id in "${fastq_files[@]}"; do
  fastq2bam "$id"
done
