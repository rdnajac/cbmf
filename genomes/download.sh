#!/bin/bash
#
# This script maintains the latest human and mouse genomes for the alignment pipelines.

# Fail on any error
set -euxo pipefail

# TODO: /mnt/data/genomes ?
LOCAL_GENOMES="$HOME/genomes"

# Create directories using brace expansion
mkdir -pv "$LOCAL_GENOMES/{human,mouse}"

# Constants for URLs and directories
NCBI="ftp.ncbi.nlm.nih.gov/genomes/all"
GENEBANK="GCA"
REFSEQ="GCF"
HUMAN_LATEST="GCA_000001405.15_GRCh38"
MOUSE_LATEST="GCA_000001635.9_GRCm39"

# Construct HTTPS URLs
HUMAN_URL="https://${NCBI}/${GENEBANK}/${HUMAN_LATEST}"
MOUSE_URL="https://${NCBI}/${GENEBANK}/${MOUSE_LATEST}"

# The folder to download the prebuilt indexes from
PREBUILT_INDEXES="seqs_for_alignment_pipelines.ucsc_ids"

# Checksum file to verify the downloaded files
CHECKSUM_FILE="md5checksums.txt"
CHECKSUMS=("${HUMAN_URL}/${CHECKSUM_FILE}" "${MOUSE_URL}/${CHECKSUM_FILE}")

# Define the flags for wget so we can easily reuse or modify them
WGET_FLAGS="--quiet --show-progress --progress=bar:force:noscroll --no-parent --no-host-directories -e robots=off --cut-dirs=7"

# Download prebuilt indexes and checksums
wget ${WGET_FLAGS} -r -P "$HOME/genomes/human" "${HUMAN_URL}/${PREBUILT_INDEXES}" "${HUMAN_URL}/${CHECKSUM_FILE}" &
wget ${WGET_FLAGS} -r -P "$HOME/genomes/mouse" "${MOUSE_URL}/${PREBUILT_INDEXES}" "${MOUSE_URL}/${CHECKSUM_FILE}" &

## Download individual checksum files
#download_checksums() {
#  for i in "${!CHECKSUMS[@]}"; do
#    wget ${WGET_FLAGS} -P "${DEST_DIRS[$i]}" "${CHECKSUMS[$i]}" &
# done
#}

## Verify checksums
#verify_checksums() {
#  for dir in "${DEST_DIRS[@]}"; do
#    pushd "$dir" > /dev/null

#    if [ -f "md5checksums.txt" ]; then
#      while read -r md5 file; do
#        if ! echo "${md5} ${file}" | md5sum -c --status; then
#          echo "Checksum verification failed for ${file}" >&2
#          exit 1
#        fi
#      done < "md5checksums.txt"
#    else
#      echo "Checksum file not found in ${dir}" >&2
#      exit 1
#    fi

#    popd > /dev/null
#  done
#}
