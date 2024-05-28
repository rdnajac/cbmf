#!/bin/bash
#
# Download and verify the latest major releases of human and mouse genomes
# and organize them into ~/genomes/{human, mouse}
#
# Author: [Your Name]
# Date: [Today's Date]

set -euo pipefail

# Constants for URLs and directories
readonly HUMAN_BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38"
readonly MOUSE_BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39"
readonly HUMAN_DIR="${HOME}/genomes/human"
readonly MOUSE_DIR="${HOME}/genomes/mouse"
readonly HUMAN_CHECKSUM_FILE="md5checksums.txt"
readonly MOUSE_CHECKSUM_FILE="md5checksums.txt"

# Helper functions for colored messages
info() {
  printf "\033[94m%s\033[0m\n" "${1}" >&2
}

warn() {
  printf "\033[91m%s\033[0m\n" "${1}" >&2
}

#######################################
# Download and verify genome files.
# Globals:
#   None
# Arguments:
#   base_url: Base URL of the genome files.
#   dest_dir: Destination directory for downloaded files.
#   checksum_file: Name of the checksum file.
# Outputs:
#   Writes status messages to STDOUT
# Returns:
#   0 if successful, exits with 1 if checksum verification fails.
#######################################
download_genome() {
  local base_url="$1"
  local dest_dir="$2"
  local checksum_file="$3"

  wget -q -r -np -nH --cut-dirs=7 -P "${dest_dir}" "${base_url}/"
  wget -q -P "${dest_dir}" "${base_url}/${checksum_file}"

  pushd "${dest_dir}" > /dev/null
  if md5sum -c "${checksum_file}"; then
    info "All files in ${dest_dir} verified successfully."
  else
    warn "Checksum verification failed for files in ${dest_dir}"
    exit 1
  fi
  popd > /dev/null
}

#######################################
# Main function to orchestrate downloads and verifications.
# Globals:
#   HUMAN_DIR
#   MOUSE_DIR
# Outputs:
#   Writes status messages to STDOUT
# Returns:
#   0 if successful
#######################################
main() {
  mkdir -p "${HUMAN_DIR}" "${MOUSE_DIR}"

  info "Downloading and verifying human genome..."
  download_genome "${HUMAN_BASE_URL}" "${HUMAN_DIR}" "${HUMAN_CHECKSUM_FILE}"

  info "Downloading and verifying mouse genome..."
  download_genome "${MOUSE_BASE_URL}" "${MOUSE_DIR}" "${MOUSE_CHECKSUM_FILE}"

  info "All downloads and verifications completed successfully."
}

main "$@"

