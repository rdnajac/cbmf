#!/bin/bash
#================================================================
## Description: Script to run FastQC on all .fastq.gz files in a
#               specified directory. It sets up an output directory
#               within the specified working directory, processes
#               each file in parallel, and waits for tasks to finish.
## Author: Ryan D. Najac
## Last modified: $(date +"%Y-%m-%d")
#================================================================

# enforce strict error handling and print each command
set -euxo pipefail

function multifastqc {
  local working_dir="$1"
  # default to 'fastq.gz' if no extension is provided
  local file_extension="${2:-fastq.gz}"
  local files=("${working_dir}"/*.$file_extension)
  local output_dir=""
  local zipfile=""

  if [[ $file_extension == "bam" ]]; then
    output_dir="${working_dir}/bamqc"
    zipfile="bamqc_html.zip"

  else
    output_dir="${working_dir}/fastqc"
    zipfile="fastqc_html.zip"
  fi

    # create output directory if it doesn't exist
    mkdir -p "$output_dir" || return 1

    # check if any files were found
    [[ ${#files[@]} -eq 0 ]] && echo "no .$file_extension files found" && return 2

    # run FastQC on each file in parallel
    fastqc -o "$output_dir" --noextract --memory 1024 -t 32 "${files[@]}"

    # compress all html files into a single zip archive
    (cd "$output_dir" && zip -r ../$zipfile *.html && rm -f *.html)
}

multifastqc "$@"

