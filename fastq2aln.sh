#!/bin/bash
#================================================================
## Description: Fastq to aligned and BAM or CRAM alignment pipeline
## Author: Ryan D. Najac
## Last modified: 2024-05-03
#================================================================
set -euo pipefail   # enforce strict error handling

readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly MAX_THREADS=$(nproc)

readonly FNA="download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz"
readonly MOUSEREF="download/GCA_000001635.9_GRCm39_genomic.fna.gz"
readonly HUMANREF="download/GRCh38_latest_genomic.fna.gz"

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
    local suffix="_001.fastq.gz"
    local r1suffix="_R1${suffix}"
    local r2suffix="_R2${suffix}"
    declare -a files=()
    local dir=${1:-$(pwd)}

    for r1 in "${dir}/"*"${r1suffix}"; do
        local sample
        sample=$(basename "$r1" | sed -e "s/${r1suffix}//")
        if [[ -f "${dir}/${sample}${r2suffix}" ]]; then
            files+=("$sample")
        fi
    done
    printf "%s\n" "${files[@]}"
}

fastq2bam() {
    local id="$1"
    local r1="${id}_R1_001.fastq.gz"
    local r2="${id}_R2_001.fastq.gz"

    bowtie2 --time --threads "$MAX_THREADS" --mm -x "$REF" -1 "$r1" -2 "$r2" \
        2>> "${id}.log" | samtools sort -@ "$MAX_THREADS" 2>> "${id}.log" > "${id}.bam"
    okay "Alignment completed for $id"
}

#======================= COMPRESSION =============================
convert() {
    local file="$1"
    case "${file##*.}" in
        bam)
            samtools view -@ "$MAX_THREADS" -C -T "$FNA" -o "${file%.bam}.cram" "$file"
            okay "Converted $file to CRAM format"
            ;;
        cram)
            samtools view -@ "$MAX_THREADS" -b -o "${file%.cram}.bam" "$file"
            okay "Converted $file to BAM format"
            ;;
        *)
            warn "Error: Unsupported file type."
            exit 1
            ;;
    esac
}

#======================= MAIN ====================================
# Example usage:
# declare -a RA=$(get_fastq_pairs ra/)
# main() {
#   for sample in ${RA[@]}; do
#     fastq2bam "${sample}"
#   done
#

