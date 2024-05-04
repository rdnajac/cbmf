#!/bin/bash
#================================================================
## Description: Additional script to generate summaries and perform
##              additional cleaning after running FastQC. It handles
##              summary report generation and file cleaning processes.
## Author: Ryan D. Najac
## Last modified: $(date +"%Y-%m-%d")
#================================================================

# Enforce strict error handling and print each command
set -euxo pipefail

# Get the script's filename for usage and debugging
SCRIPT_NAME=$(basename "$BASH_SOURCE" .sh)

# Exit function with error message and status code
bail() {
    echo -ne "$1" >&2
    exit ${2:-1}
}

HELP_MSG="Usage: $SCRIPT_NAME <working_directory>\n
Options:
  -h    Display this help and exit
Example:
  $SCRIPT_NAME /path/to/data   # Run summary and cleaning scripts in /path/to/data
"

usage() {
    local status=2
    if [ "$1" -eq "$1" ] 2>/dev/null; then
        status=$1
        shift
    fi
    bail "${1}${HELP_MSG}" $status
}

# Parse command-line options
while getopts "h" opt; do
    case $opt in
        h) usage 0 ;;
        ?) usage "Invalid option: -$OPTARG \n" ;;
    esac
done

# Validate the number of arguments
shift $((OPTIND - 1))
if [[ $# -lt 1 ]]; then
    usage "Too few arguments\n"
fi

# Set the path to the qcsummary.py to the same directory as this script
QCSUMMARY_PY=$(dirname "$0")/qcsummary.py

# Ensure the qcsummary.py script exists and is executable
[[ -f $QCSUMMARY_PY ]] || bail "Error: $QCSUMMARY_PY not found\n"
[[ -x $QCSUMMARY_PY ]] || chmod +x $QCSUMMARY_PY

# Run the qcsummary.py script to generate a summary report
./$QCSUMMARY_PY "$1"

# Ensure the performance_metrics.csv file exists
PERFORMANCE_METRICS_CSV="performance_metrics.csv"
[[ -f $PERFORMANCE_METRICS_CSV ]] || bail "Error: $PERFORMANCE_METRICS_CSV not found\n"

ADDITIONAL_CLEANING_PY="$(dirname "$0")/additional_cleaning.py"
[[ -f $ADDITIONAL_CLEANING_PY ]] || bail "Error: $ADDITIONAL_CLEANING_PY not found\n"
[[ -x $ADDITIONAL_CLEANING_PY ]] || chmod +x $ADDITIONAL_CLEANING_PY

# Run the additional_cleaning.py script
$ADDITIONAL_CLEANING_PY "$PERFORMANCE_METRICS_CSV" "./flagstat.csv"


