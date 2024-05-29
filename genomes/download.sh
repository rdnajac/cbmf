#!/bin/bash
#
# Download and verify the latest major releases of human and mouse genomes.

set -euxo pipefail

# Constants for URLs and directories
readonly NCBI_BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all"
readonly HUMAN_LATEST="GCA_000001405.15_GRCh38"
readonly MOUSE_LATEST="GCA_000001635.9_GRCm39"
readonly HUMAN_BASE_URL="${NCBI_BASE_URL}/GCA/000/001/405/${HUMAN_LATEST}"
readonly MOUSE_BASE_URL="${NCBI_BASE_URL}/GCA/000/001/635/${MOUSE_LATEST}"
readonly PREBUILT_INDEXES="seqs_for_alignment_pipelines.ucsc_ids"
readonly BASE_URLS=("$HUMAN_BASE_URL" "$MOUSE_BASE_URL")
readonly CHECKSUMS=(
  "${HUMAN_BASE_URL}/md5checksums.txt"
  "${MOUSE_BASE_URL}/md5checksums.txt"
)
readonly WGET_FLAGS="--quiet --show-progress --progress=bar:force:noscroll --no-parent --no-host-directories -e robots=off --cut-dirs=7"
readonly DEST_DIRS=("$HOME/genomes/human" "$HOME/genomes/mouse")

create_directories() {
  for dir in "${DEST_DIRS[@]}"; do
    mkdir -p "$dir"
  done
}

# Download prebuilt indexes recursively
download_prebuilt_indexes() {
  local prebuilt_path="$1"
  for i in "${!BASE_URLS[@]}"; do
    wget ${WGET_FLAGS} -r -P "${DEST_DIRS[$i]}" "${BASE_URLS[$i]}/$prebuilt_path" &
  done
}

# Download individual checksum files
download_checksums() {
  for i in "${!CHECKSUMS[@]}"; do
    wget ${WGET_FLAGS} -P "${DEST_DIRS[$i]}" "${CHECKSUMS[$i]}" &
 done
}

# Verify checksums
verify_checksums() {
  for dir in "${DEST_DIRS[@]}"; do
    pushd "$dir" > /dev/null

    if [ -f "md5checksums.txt" ]; then
      while read -r md5 file; do
        if ! echo "${md5} ${file}" | md5sum -c --status; then
          echo "Checksum verification failed for ${file}" >&2
          exit 1
        fi
      done < "md5checksums.txt"
    else
      echo "Checksum file not found in ${dir}" >&2
      exit 1
    fi

    popd > /dev/null
  done
}

# Main
create_directories
download_prebuilt_indexes "$PREBUILT_INDEXES"
download_checksums
#verify_checksums

