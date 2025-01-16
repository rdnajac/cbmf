#!/bin/bash

#!/bin/bash

align_with_hisat2() {
	local sampleID="$1"
	local INPUT_FOLDER="$2"
	local OUTPUT_FOLDER="$3"
	local THREADS="$4"
	local REFERENCE="$5"

	# Derive the R1 and R2 FASTQ file names from sampleID
	local r1_fastq="${INPUT_FOLDER}/${sampleID}_R1*.fastq.gz"
	local r2_fastq="${INPUT_FOLDER}/${sampleID}_R2*.fastq.gz"

	if [[ ! -f "$r1_fastq" || ! -f "$r2_fastq" ]]; then
		echo "Error: Missing R1 or R2 file for sample $sampleID"
		return 1
	fi

	# Create output BAM file name
	local bam_file="${OUTPUT_FOLDER}/${sampleID}_sorted_markdup.bam"

	# Align reads and process BAM file
	# the `dta` flag is only if we're aligning RNAseq reads
	hisat2 -p "$THREADS" --dta -x "$REFERENCE" -1 "$r1_fastq" -2 "$r2_fastq" \
		| samtools sort -n -@ "$THREADS" - \
		| samtools fixmate -@ "$THREADS" -m - - \
		| samtools sort -@ "$THREADS" - \
		| samtools markdup -@ "$THREADS" - "$bam_file"

	echo "Alignment completed for $sampleID. Output BAM: $bam_file"
}

# Example usage:
# align_with_hisat2 "sampleID" "/path/to/input" "/path/to/output" 4 "/path/to/reference"
