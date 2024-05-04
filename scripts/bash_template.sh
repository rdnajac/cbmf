#================================================================
## Description: fastq to aligned and bam or cram alignment pipeline
## Author: Ryan D. Najac
## Last modified: 2024-05-03
#================================================================

#!/bin/bash

#======================= SET UP =================================

# TODO toggle this with a flag
set -euo pipefai        # enforce strict error handling
set -x                   # print each command before execution

# over-engineered way to get the script directory...
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCRIPT_NAME=$(basename "$0")

# exit and print error message with status code
bail() { echo -ne "$1" >&2; exit ${2:-1}; }


#======================= HELP ===================================

HELP_MSG="Usage: $SCRIPT_NAME <directory> [file_extension]\n
Options:
  -h    Display this help and exit
Example:
  $SCRIPT_NAME /path/to/data            # Process all files in /path/to/data
"

usage() {
    local status=2
    [ "$1" -eq "$1" ] 2>/dev/null && status=$1 && shift
    bail "${1}${HELP_MSG}" $status
}


#======================= COMMAND LINE PARSING ===================

while getopts "h" opt; do
    case $opt in
        h) usage 0 ;;
        ?) usage "Invalid option: -$OPTARG \n" ;;
    esac
done

shift $((OPTIND - 1)) && [[ $# -lt 1 ]] && usage "Too few arguments\n"

#==========MAIN CODE BELOW==========

