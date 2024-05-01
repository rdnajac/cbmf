#!/bin/bash
#================================================================
## Description: Script to run FastQC on all .fastq.gz files in a
#               specified directory. It sets up an output directory
#               within the specified working directory, processes
#               each file in parallel, and waits for tasks to finish.
## Author: Ryan D. Najac
#================================================================

# strict enforcement of error handling
set -euxo pipefail

# get the script name from filename
SCRIPT_NAME=$(basename "$BASH_SOURCE" .sh)

# exit with an error message and status code
bail() { echo -ne "$1" >&2; exit ${2:-1}; }

# usage message and exit
HELP_MSG="Usage: $SCRIPT_NAME <directory> [file_extension]
Options:
  -h    display this help and exit
"

usage() {
    local status=2
    if [ "$1" -eq "$1" ] 2>/dev/null; then
        status=$1
        shift
    fi
    bail "${1}${HELP_MSG}" $status
}

# option parsing
while getopts "h" opt; do
    case $opt in
        h) usage 0 ;;
        ?) usage "Invalid option: -$OPTARG \n" ;;
    esac
done

# check for minimum required arguments
shift $((OPTIND - 1))
if [[ $# -lt 1 ]]; then
    usage "Too few arguments\n"
fi

#==========MAIN CODE BELOW==========

function multifastqc {
    local working_dir="$1"
    # default to 'fastq.gz' if no extension is provided
    local file_extension="${2:-fastq.gz}"
    local output_dir=""
    if [[ $file_extension == "bam" ]]; then
        output_dir="${working_dir}/bamqc"
    else
        output_dir="${working_dir}/fastqc"
    fi
    local files=("${working_dir}"/*.$file_extension)

    # create output directory if it doesn't exist
    mkdir -p "$output_dir" || return 1

    # check if any files were found
    [[ ${#files[@]} -eq 0 ]] && echo "no .$file_extension files found" && return 2

    # run FastQC on each file in parallel
    for file in "${files[@]}"; do
        fastqc -o "$output_dir" --noextract --memory 1024 "$file" &
    done
    wait # for all background jobs to finish
}

# run the multifastqc function with the provided arguments
multifastqc "$@"
