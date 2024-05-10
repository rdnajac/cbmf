#!/bin/bash
#================================================================
## Description: Fastq to aligned and BAM or CRAM alignment pipeline
## Author: Ryan D. Najac
## Last modified: 2024-05-03
#================================================================
set -euo pipefail   # enforce strict error handling
set -x

readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly THREADS=$(nproc)

readonly MOUSEREF="download/GCA_000001635.9_GRCm39_genomic.fna.gz"
readonly HUMANREF="/home/ubuntu/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"

#======================= HELP AND OPTIONS =======================
readonly OPTIONS="hxg:"
readonly USAGE="Usage: ${SCRIPT_NAME} [-hxg] <directory>
Options:
-h    Display this help and exit
-x    Print each command before execution
-g    Specify reference genome ('mouse' or 'human', default: 'human')
"

# Helper functions to print colored messages
okay() { printf "\033[92m%s\033[0m\n" "${1}" >&2; }
info() { printf "\033[94m%s\033[0m\n" "${1}" >&2; }
warn() { printf "\033[91m%s\033[0m\n" "${1}" >&2; }

# Quit with an error message
bail () { warn "${1}"; exit ${2:-1}; }

# Print usage and exit
usage(){ printf "%s\n"; "$USAGE"; exit ${1:-0} >&2; }

#======================= COMMAND LINE PARSING ===================
declare REF_GENOME="human"  # Default to human genome
while getopts ":hxg:" opt; do
  case $opt in
    h) usage 0 ;;
    x) set -x ;;
    g) REF_GENOME=${OPTARG} ;;
    \?) printf "Invalid option: -%s\n" "$OPTARG" >&2; usage 1 ;;
  esac
done

shift $((OPTIND - 1))
if [[ $# -lt 1 ]]; then
  warn "Error: Too few arguments"
  usage 1
fi

# Select reference genome based on the option
case "$REF_GENOME" in
  mouse) REF="$MOUSEREF" ;;
  human) REF="$HUMANREF" ;;
  *)     warn "Error: Unsupported reference genome '$REF_GENOME'."; exit 1 ;;
esac

#======================= MAIN FUNCTIONALITY =====================
get_fastq_pairs() {
  local dir=${1:-$(pwd)}
  local suffix="_001.fastq.gz"
  declare -a files=()

    # Use find to avoid globbing issues and directly process each matching file
    while IFS= read -r r1; do
      local r2="${r1%_R1$suffix}_R2$suffix"
      if [[ -f "$r2" ]]; then
        files+=("$(basename "${r1%_R1$suffix}")")
      fi
    done < <(find "$dir" -maxdepth 1 -type f -name "*_R1$suffix")

    printf "%s\n" "${files[@]}"
  }

fastq2bam() {
    local id="$1"
    local r1="${id}_R1_001.fastq.gz"
    local r2="${id}_R2_001.fastq.gz"

    # Align the reads, sort by name, fix mate information, sort again, and mark duplicates
    bowtie2 --time --threads "$THREADS" --mm -x "$REF" -1 "$r1" -2 "$r2" 2>> "${id}.log" \
            | samtools sort -n -@ "$THREADS" 2>> "${id}.log" > "${id}.bam"
      # | samtools sort -n -@ "$THREADS" - \
      # | samtools fixmate -m -@ "$THREADS" - - \
      # | samtools sort -@ "$THREADS" - \
      # | samtools markdup -r -@ "$THREADS" - "${id}_dedup.bam"
}

#======================= MAIN ====================================
#[[ s -d "$1" ]] && bail "Error: Directory '$1' not found."

declare -a fastq_files=()
fastq_files=($(get_fastq_pairs "$1"))

for id in "${fastq_files[@]}"; do
  fastq2bam "$id"
done
