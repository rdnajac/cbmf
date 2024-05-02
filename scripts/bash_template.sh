#!/bin/bash
#================================================================
## Description: This script wraps other scripts
## Author: Ryan D. Najac
## Last modified: $(date +"%Y-%m-%d")
#================================================================

# enforce strict error handling and print each command
set -euxo pipefail

# get the script name from filename
SCRIPT_NAME=$(basename "$BASH_SOURCE" .sh)

# exit and print error message with status code
bail() { echo -ne "$1" >&2; exit ${2:-1}; }

HELP_MSG="Usage: $SCRIPT_NAME <directory> [file_extension]\n
Options:
  -h    Display this help and exit
Example:
  $SCRIPT_NAME /path/to/data            # Process all files in /path/to/data
"

usage() {
    local status=2
    if [ "$1" -eq "$1" ] 2>/dev/null; then
        status=$1
        shift
    fi
    bail "${1}${HELP_MSG}" $status
}

# parse command-line options
while getopts "h" opt; do
    case $opt in
        h) usage 0 ;;
        ?) usage "Invalid option: -$OPTARG \n" ;;
    esac
done

# validate argument count
shift $((OPTIND - 1))
if [[ $# -lt 1 ]]; then
    usage "Too few arguments\n"
fi

#==========MAIN CODE BELOW==========

