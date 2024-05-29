#!/bin/bash
#================================================================
## Description: Template for bash scripts
## Author: Ryan D. Najac
## Last Modified: 2024-05-03
#================================================================

#======================= SET UP ==================================
set -euo pipefail   # Enforce strict error handling
set -x              # Print each command before execution

# Secure and robust way to get the script directory
readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")

# Exit and print error message with status code
bail() { printf '%s\n' "$1" >&2; exit "${2:-1}"; }

#======================= HELP ====================================
readonly HELP_MSG="
Usage: $SCRIPT_NAME <directory> [file_extension]
Options:
  -h    Display this help and exit
Example:
  $SCRIPT_NAME /path/to/data            # Process all files in /path/to/data
"

# Print help message and exit
usage() { local status=${1:-2}; printf '%s\n' "$HELP_MSG"; exit "$status"; }

#======================= COMMAND LINE PARSING ====================
while getopts ":h" opt; do
    case "$opt" in
        h) usage 0 ;;
        \?) usage "Invalid option: -$OPTARG\n" 1 ;;
    esac
done

shift $((OPTIND - 1))

# Check if at least one argument is provided
if [[ $# -lt 1 ]]; then
    usage "Too few arguments\n" 1
fi

DIRECTORY=$1
FILE_EXTENSION=${2:-*}  # Default to all file types if not specified

#========== MAIN CODE BELOW ==========
printf "Processing files in %s with extension %s\n" "$DIRECTORY" "$FILE_EXTENSION"
# Example command to demonstrate processing
find "$DIRECTORY" -type f -name "*.$FILE_EXTENSION" -exec printf "Processing file: %s\n" {} \;

