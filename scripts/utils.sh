#!/bin/bash
#
# Source this file to configure the script environment
# https://gabrielstaples.com/bash-libraries
# https://superuser.com/questions/46139/what-does-source-do
# https://askubuntu.com/questions/182012/is-there-a-difference-between-and-source-in-bash-after-all

# Check if the script has already been sourced
if [ -z "${UTILS_SOURCED}" ]; then
    UTILS_SOURCED=1
else
    return
fi

export GENOMES="/home/ubuntu/genomes"
export MOUSEREF="${GENOMES}/mouse/GCA_000001635.9_GRCm39_full_analysis_set"
export HISAT_MOUSE_INDEX="${MOUSEREF}.fna.hisat2_index"
export MREF_GTF="${MOUSEREF}.refseq_annotation.gtf"
export MREF_GFF="${MOUSEREF}.refseq_annotation.gff"

# Define ANSI color codes for colored output
export RED='\033[0;31m'
export GRN='\033[0;32m'
export YEL='\033[0;33m'
export BLU='\033[0;34m' 
export MAG='\033[0;35m'
export CYN='\033[0;36m'
export WHT='\033[0;37m'
export RESET='\033[0m'

# Helper functions to print colored messages
warn() { printf "${RED}%s${RESET}\n" "$1" >&2; }
info() { printf "${BLU}%s${RESET}\n" "$1"; }

# Quit with an error message and optional exit status
bail() { warn "$1"; exit "${2:-1}"; }

# Error Codes
export E_COMMAND_NOT_EXECUTABLE=126
export E_COMMAND_NOT_FOUND=127

# Check if a command is available
ensure_installed() {
    if ! command -v "$1" &> /dev/null; then
        bail "Error: $1 is not installed or not in PATH" "$E_COMMAND_NOT_FOUND"
    fi
}

# If this script is being executed, warn the user
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    warn "This script is not meant to be executed directly."
    exit 1
else 
    info "Sourced ${BASH_SOURCE[0]}"
fi
