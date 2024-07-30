#!/bin/bash
#
## Demultplex Illumina Sequencing Data with bcl2fastq
## note: this script assumes that the run folder contains
## a valid SampleSheet.csv
## all par

demultiplex_bcl2fastq()
{
	local run_folder=$1
	local output_folder=$2

	# run bcl2fastq with the same flags as a default NextSeq run
	bcl2fastq --no-lane-splitting \
		--ignore-missing-bcls \
		--ignore-missing-filter \
		--ignore-missing-positions \
		--ignore-missing-controls \
		--auto-set-to-zero-barcode-mismatches \
		--find-adapters-with-sliding-window \
		--adapter-stringency 0.9 \
		--mask-short-adapter-reads 35 \
		--minimum-trimmed-read-length 35 \
		-R "$run_folder" \
		-o "$output_folder"
	# --sample-sheet "$sample_sheet" \

	echo "Demultiplexing complete"
	exit 0
}

# check if the correct number of arguments were provided
if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <run_folder> <output_folder>"
	exit 1
fi

demultiplex_bcl2fastq "$1" "$2"
