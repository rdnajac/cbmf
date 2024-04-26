#!/bin/bash

#***********************************************************************#
#                           multi_fastqc.sh                             #
#                   written by [Your Name Here]                         #
#                      [Date Last Modified]                             #
#                                                                       #
#       This script runs FastQC on all FASTQ files in the current       #
#       directory and outputs the results to a specified directory.     #
#***********************************************************************#

# Exit codes
E_BADDIR=85  # No such directory

# Run FastQC on multiple FASTQ files
# Globals:
#   None
# Arguments:
#   None
# Outputs:
#   Creates files in the fastqc_results directory
# Returns:
#   E_BADDIR if output directory cannot be created
multi_fastqc() {
    local output_dir="./fastqc_results"
    mkdir -p "$output_dir" || return $E_BADDIR

    local fastqs=(./*.fastq.gz)
    echo "Running FastQC on files in the current directory..."

    # Iterate over each FASTQ file and process it with FastQC
    for fastq in "${fastqs[@]}"; do
        echo "Processing $fastq..."
        fastqc -o "$output_dir" --noextract --memory 1024 "$fastq" &
    done

    wait # Wait for all background jobs to finish
    echo "FastQC processing complete."
}

# Check if script is being run directly and if so, call multi_fastqc
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    multi_fastqc
fi

exit 0

