#!/bin/bash
## Script to convert BCL files to FASTQ format using bcl2fastq
## Note that this script assumes there is a valid `samplesheet.csv` in the run folder

# Exit immediately if a command exits with a non-zero status
set -e

# Function to display usage information
usage() {
	echo "Usage: $0 -r <run_folder> -o <output_folder>"
	echo
	echo "Options:"
	echo "  -r  Path to the run folder containing BCL files"
	echo "  -o  Path to the output folder for FASTQ files"
	exit 1
}

# Parse command line arguments
while getopts ":r:o:" opt; do
	case $opt in
	r) run_folder="$OPTARG" ;;
	o) output_folder="$OPTARG" ;;
	\?)
		echo "Invalid option -$OPTARG" >&2
		usage
		;;
	esac
done

# Check if required arguments are provided
if [ "$run_folder" = "" ] || [ "$output_folder" = "" ]; then
	echo "Error: Both run folder and output folder must be specified."
	usage
fi

# Check if the run folder exists
if [ ! -d "$run_folder" ]; then
	echo "Error: Run folder does not exist: $run_folder"
	exit 1
fi

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Run bcl2fastq command
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

echo "BCL to FASTQ conversion completed successfully."
