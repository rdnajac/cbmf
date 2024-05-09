#!/bin/bash
#================================================================
## Description: Script to run FastQC on various types of sequence files in a
## specified input directory. Outputs are generated in a specified output directory.
## Author: Ryan D. Najac
## Date of modification: $(date +"%Y-%m-%d")
#================================================================
set -euo pipefail

function multifastqc {
    local input_dir="$1"
    local output_dir="$2"

    # Create output directory if it doesn't exist
    mkdir -pv "$output_dir" || return 1

    # Explicit remote-compatible command to find files
    local files=( $(ssh aws "ls ${input_dir}/*.fastq ${input_dir}/*.fastq.gz ${input_dir}/*.bam 2>/dev/null") )

    # Check if any files were found
    if [ ${#files[@]} -eq 0 ]; then
        echo "No .fastq, .fastq.gz, or .bam files found in ${input_dir}"
        return 2
    fi

    # Determine the appropriate file extension for naming the zip file
    local extension="${files[0]##*.}"
    local zipfile="${output_dir}/${extension}qc_html.zip"

    # Run FastQC on each file in parallel
    ssh aws "fastqc -o ${output_dir} --noextract --memory 1024 -t \$(nproc) ${files[*]}"

    # Zip HTML results and clean up
    ssh aws "zip -r ${output_dir}/${zipfile} ${output_dir}/*.html && rm -f ${output_dir}/*.html"
    # this didnt work
}

readonly PATH_TO_FILES="/home/ubuntu/rnaseq/cloneB"
multifastqc  "$PATH_TO_FILES/fastq" "fastqc"
multifastqc  "$PATH_TO_FILES/aligned" "bamqc"

