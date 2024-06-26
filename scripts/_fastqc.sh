#!/bin/bash

# Arbitrarily setting the memory limit to 4GB
MEMSIZE=4096

#######################################
# Batch processing of FastQC
# Get configuration directory.
# Globals:
#   MEMSIZE
# Arguments:
#  $1: Input directory
#  $2: Output directory (optional; Default: $1/fastqc)
# Outputs:
#   Writes location to stdout
#
# Details:
# fastqc attemts to run on every file in the input directory and,
# by default, its subdirectories. It will skip any file with an
# extension that is not recognized by fastqc (ie .fastq.gz, .bam)
#######################################
_fastqc() {
	local input_dir="$1"
	local output_dir="${2:-fastqc}"

	cd "$input_dir" || {
		echo "Error: $input_dir not found" >&2
		exit 1
	}

	# Create output directory
	mkdir -vp "$output_dir"

	## Run FastQC
	# --noextract: Do not extract the results
	# -t: Number of threads
	# --memory: Memory limit for java processes
	# fastqc --noextract -t "$(nproc)" --memory "$MEMSIZE" -o "$output_dir" ./*
	fastqc --noextract -t "$(nproc)" --memory 4096 -o bamqc ./*

}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && _fastqc "$@"
