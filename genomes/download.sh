#!/bin/bash
#
# This script maintains the latest human and mouse genomes for the alignment pipelines.

# Fail on any error
set -euxo pipefail

# Base URL and constants for URLs and directories
NCBI_FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all"
HUMAN_GENOME_PATH="GCA/000/001/405/GCA_000001405.15_GRCh38"
MOUSE_GENOME_PATH="GCA/000/001/635/GCA_000001635.9_GRCm39"
PREBUILT_INDEXES="seqs_for_alignment_pipelines.ucsc_ids"

# Construct full URLs
HUMAN_URL="${NCBI_FTP_BASE}/${HUMAN_GENOME_PATH}"
MOUSE_URL="${NCBI_FTP_BASE}/${MOUSE_GENOME_PATH}"

# Define the flags for wget so we can easily reuse or modify them
WGET_NO_FLAGS="--no-parent --no-directories --no-host-directories --no-check-certificate"
WGET_ADDITIONAL_FLAGS="--quiet -e robots=off --cut-dirs=7 --show-progress --progress=bar:force:noscroll"
MY_WGET="wget ${WGET_NO_FLAGS} ${WGET_ADDITIONAL_FLAGS}"

download_genome() {
  local url="$1"
  local species="$2"

  mkdir -pv "$species"

  # Download the genome files in parallel
  for file in $(wget "${url}" -O - | grep -oP '(?<=href=")[^"]*' | grep -E '(\.bed|\.fai|\.gz|\.gz\.txt)$'); do
    ${MY_WGET} -P "$species" "${url}/${file}" &
  done
}

if [[ "$1" =~ ^(-h|--human|[Hh][Uu][Mm][Aa][Nn])$ ]]; then
  download_genome "${HUMAN_URL}/${PREBUILT_INDEXES}" "human"
elif [[ "$1" =~ ^(-m|--mouse|[Mm][Oo][Uu][Ss][Ee])$ ]]; then
  download_genome "${MOUSE_URL}/${PREBUILT_INDEXES}" "mouse"
else
  echo "Invalid argument: $1" >&2
  exit 1
fi

echo -e "🌈\e[31mS\e[32mu\e[33mc\e[34mc\e[35me\e[36ms\e[37ms\e[0m!✨"

# TODO check checksums
