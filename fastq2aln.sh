#================================================================
## Description: fastq to aligned and bam or cram alignment pipeline
## Author: Ryan D. Najac
## Last modified: 2024-05-03
#================================================================

#!/bin/bash

#======================= SET UP =================================
SCRIPT_NAME=$(basename "$0")
SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MAX_THREADS=$(nproc)

set -euo pipefail   # enforce strict error handling
# set -x              # print each command before execution
                    # TODO toggle these options with a flag

# exit and print error message with status code
bail() { echo -ne "$1" >&2; exit ${2:-1}; }

#======================= DEFAULTS ===============================
FNA="download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz"
# MOUSEREF="download/GCA_000001635.9_GRCm39_genomic.fna.gz"
# HUMANREF="download/GRCh38_latest_genomic.fna.gz"
#
#======================= HELP ===================================
readonly OPTIONS=":hxm"
# TODO long options?
readonly USAGE="Usage: ${SCRIPT_NAME} [${OPTIONS}] <directory>
options:
  -h    display this help and exit
  -x    print each command before execution

  --mouse    use the mouse reference genome
  --human    use the human reference genome
"

usage_and_exit()
{
    local status=2
    [ "$1" -eq "$1" ] 2>/dev/null && status=$1 && shift
    bail "${1}${USAGE}" $status
}

#======================= COMMAND LINE PARSING ===================
# see `man getopts`
while getopts "h" opt; do
    case $opt in
        h) usage_and_exit 0 ;;
        ?) usage_and_exit "Invalid option: -$OPTARG \n" ;;
    esac
done
# make sure the script has at least one argument
# TODO is there a better way to do this?
shift $((OPTIND - 1)) && [[ $# -lt 1 ]] && usage "Too few arguments\n"

#======================= MD5 ====================================

# Generate checksums for all files in the current directory
generate_checksums()
{
    local extension=${1:-"fastq.gz"}
    local checksum_file="md5sums.txt"
    echo "Generating checksums for *.$extension files in the current directory..."
    for file in *.$extension; do
        echo "$(md5sum "$file" | cut -d ' ' -f 1)\t$file" >> $checksum_file
    done
    echo "Checksums stored in $checksum_file"
}

# Verify checksums for all files in the current directory
verify_checksums()
{
    local md5sums_txt="${1:-md5sums.txt}"
    echo "Verifying checksums..."
    while read -r expected_checksum file; do
        local status="[ OK ]"
        [[ "$(md5sum "$file" | cut -d ' ' -f 1)" != "$expected_checksum" ]] && status="[FAIL]"
        echo "${status} ${file}"
    done < $md5sums_txt
}

#======================= UTILITY ===============================

# Function: get_fastq_pairs
# Usage: get_fastq_pairs <directory>
#
# Description:
# get a newline-delimited list of fastq R1/R2 pairs as
# sanitized IDs that correspond to the original sample ID
get_fastq_pairs ()
{
  local suffix="_001.fastq.gz"
  local r1suffix="_R1${suffix}"
  local r2suffix="_R2${suffix}"
  declare -a files=()
  local dir=${1:-$(pwd)}
  [[ "${dir}" == */ ]] && dir="${dir%/}"  # remove trailing `/`

  for r1 in "${dir}/"*"${r1suffix}"; do
    # check if we have a matching pair and add to array
    sample=$(basename "$r1" | sed -e "s/${r1suffix}//")
    [[ -f "${dir}/${sample}${r2suffix}" ]]; && files+=("$sample")
  done
  printf "%s\n" "${files[@]}"
}

#======================= ALIGNMENT ==============================
readonly REF="$MOUSEREF"  # TODO make this an option

# Run bowtie2 with the specified reference genome using the
# specified number of threads and log the time taken
readonly BT2="bowtie2 --time --threads ${THREADS} --mm -x ${REF}"

##
# Function: fastq2bam
# Description: Align paired-end reads to the reference genome
#
# Details:
# Align reads with bowtie2 and convert to BAM format
# Declaring this as a function lets us use it in a loop
fastq2bam()
{
  local r1="${1}_R1_001.fastq.gz" && [[ ! -f "${r1}" ]] && bail "error: ${r1} not found" 2
  local r2="${1}_R2_001.fastq.gz" && [[ ! -f "${r2}" ]] && echo "error: ${r2} not found" 2
  # TODO do I need the bail function and all this error checking if I'm using `set -e`?
  $BT2 -1 "${r1}" -2 "${r2}" 2>> "${id}.log" | samtools sort -@ "${THREADS}" 2>> "${id}.log" > "${id}.bam"
}


#======================= COMPRESSION =============================
# CRAM is a compressed BAM format that is more efficient.

readonly CONVERT="samtools view -@ $THREADS"

# Helper functions with automatic filetype detection
embiggen () { ${CONVERT} -b         -o ${1%.cram}.bam $1; }
unbiggen () { ${CONVERT} -C -T $FNA -o ${1%.bam}.cram $1; }

# Function: convert
# Usage: convert <file>
# Description: Convert a BAM file to CRAM or vice versa
convert()
{
  case "${1##*.}" in
    bam) unbiggen "$1" ;;
    cram) embiggen "$1" ;;
    *) echo "error: Unsupported file type."; exit 1 ;;
  esac
}

#======================= MAIN ==================================

#================================================================
# Example usage:
# declare -a RA=$(get_fastq_pairs ra/)
# main() {
#   for sample in ${RA[@]}; do
#     fastq2bam "${REFERENCE_GENOME}" "${sample}"
#   done
# }
# main

