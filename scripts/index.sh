#!/bin/bash
#================================================================
## get and manipulate a list of fastq files
## useful for processing multiple fastq files in a directory
#================================================================

function usage() {
    echo "Usage: $0 [<directory>]"
    exit 0
}

function get_fastq_pairs {
    local dir=${1:-$(pwd)}
    declare -a files=()
    for r1 in "${dir}"/*_R1_001.fastq.gz; do
        local r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
        [[ -f "$r2" ]] && files+=("${r1%_R1_001.fastq.gz}") || \
          echo "Error: R2 file not found for ${r1##*/}"
    done
    printf "%s\n" "${files[@]}"
}

function test_get_fastq_files {
    local test_dir="tmp"
    echo "Creating temporary directory and test files..."
    mkdir -p "$test_dir"
    touch "$test_dir"/{A..Z}_R{1..2}_001.fastq.gz

    echo "Testing get_fastq_files function with simulated fastq files:"
    get_fastq_files "$test_dir"

    echo "Cleaning up test files and directory..."
    rm -r "$test_dir"
}

get_fastq_pairs $1
