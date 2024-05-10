#!/bin/bash
#================================================================
## Description: Script to run FastQC on various types of sequence files in a
## specified input directory. Outputs are generated in a specified output directory.
## Author: Ryan D. Najac
## Date of modification: $(date +"%Y-%m-%d")
#================================================================
set -euo pipefail
set -x
function multifastqc {
    local input_dir="$1"
    local output_dir="$2"
    local zipped_html_file="$3"

    mkdir -pv "$output_dir" || return 1

    fastqc -o "$output_dir" --noextract --memory 1024 -t "$(nproc)" "$input_dir"/*

    [[ -n "$zipped_html_file" ]] && zip -r "$zipped_html_file" "$output_dir"/*.html
}

readonly PATH_TO_FILES="/home/ubuntu/rnaseq/"
readonly PROJECT_FOLDERS=("agx51"          "ra")

for folder in "${PROJECT_FOLDERS[@]}"; do
    fastqc_output="${PATH_TO_FILES}/${folder}/fastqc"
    bamqc_output="${PATH_TO_FILES}/${folder}/bamqc"
    multifastqc "${PATH_TO_FILES}/${folder}/fastq" "$fastqc_output" "${fastqc_output}/fastqc_html.zip" &
    fastqc_pid=$!

    multifastqc "${PATH_TO_FILES}/${folder}/aligned" "$bamqc_output" "${bamqc_output}/bamqc_html.zip"

    wait $fastqc_pid
    echo "Finished processing $folder"
done

