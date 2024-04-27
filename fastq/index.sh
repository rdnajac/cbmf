#!/bin/bash

# Function to get a list of all `.fastq.gz` files in the specified directory
# If they're paired-end reads, they should have the same prefix, and the cleaned array is returned
function get_fastq_files {
    local dir=$1
    declare -a r1_files=()
    declare -a r2_files=()
    declare -a paired_files=()

    # Collect R1 and R2 files into separate arrays
    for file in "${dir}"/*.fastq.gz; do
        if [[ "$file" =~ _R1_001.fastq.gz ]]; then
            r1_files+=("$file")
        elif [[ "$file" =~ _R2_001.fastq.gz ]]; then
            r2_files+=("$file")
        fi
    done

    # Match R1 and R2 files based on prefix and remove suffix
    for r1 in "${r1_files[@]}"; do
        local prefix=${r1%%_R1_001.fastq.gz}
        local match=0
        for r2 in "${r2_files[@]}"; do
            if [[ "$r2" == "${prefix}_R2_001.fastq.gz" ]]; then
                paired_files+=("${prefix##*/}")
                match=1
                break
            fi
        done
        if [ "$match" -eq 0 ]; then
            echo "Error: R2 file not found for ${prefix##*/}"
            exit 1
        fi
    done

    echo "Paired prefixes: ${paired_files[@]}"
}

# Test function to create a temporary directory, simulate fastq files, and test the get_fastq_files function
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

# Run the test function
test_get_fastq_files

