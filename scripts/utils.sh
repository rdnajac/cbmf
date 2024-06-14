#!/bin/bash
#
# Source this file to configure the script environment
# https://superuser.com/questions/46139/what-does-source-do
# https://askubuntu.com/questions/182012/is-there-a-difference-between-and-source-in-bash-after-all


# Global variables
export THREADS=$(nproc)

# Path to the reference genome
export GENOMES="/home/ubuntu/genomes"
export MOUSEREF="${GENOMES}/mouse/GCA_000001635.9_GRCm39_full_analysis_set"
export HISAT_MOUSE_INDEX="${MOUSEREF}.fna.hisat2_index"
export MREF_GTF="${MOUSEREF}.refseq_annotation.gtf"
export MREF_GFF="${MOUSEREF}.refseq_annotation.gff"

# Global settings _____________________________________________________ 

set -euo pipefail   # enforce strict error handling


SCRIPT_NAME=$(basename "${0}")
USAGE="Usage: ${SCRIPT_NAME}"

# Helper functions to print colored messages
warn() { printf "\033[91m%s\033[0m\n" "${1}" >&2; }
info() { printf "\033[92m%s\033[0m\n" "${1}" >&2; }

# Quit with an error message
bail () { warn "${1}"; exit "${2:-1}"; }

# Print usage and exit
usage(){ printf "%s\n"; "$USAGE"; exit "${1:-0}" >&2; }

