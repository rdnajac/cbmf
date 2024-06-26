#!/bin/bash
#
## Wrapper script for paired-end read alignment using HISAT2 or Bowtie2.
## Calls _align() to run either bowtie2 or hisat2 on a pair of fastq files
set -o pipefail
set -x

# THREADS="$(nproc)"

_align() {

	local ALIGNER="$1"
	local R1_FASTQ="$2"
	local R2_FASTQ="$3"
	local REFERENCE="$4"

	if [ "$ALIGNER" == "hisat2" ]; then
		hisat2 -p "$(nproc)" --mm -1 "$R1_FASTQ" -2 "$R2_FASTQ" -x "$REFERENCE --dta"
		# add flag for downstream transcriptome assembly
	elif [ "$ALIGNER" == "bowtie2" ]; then
		bowtie2 -p "$(nproc)" --mm -1 "$R1_FASTQ" -2 "$R2_FASTQ" -x "$REFERENCE"
	else
		echo "Invalid aligner: $ALIGNER"
		exit 1
	fi
}

# run align and capture the output. pipe stdout to samtools and stderr to log file
_align2bam() {
	local align_cmd="$1"
	local r1_fastq="$2"
	local r2_fastq="$3"
	local reference="$4"
	local log_file="${r1_fastq%_R1*}_align.log"
	local bam_file="${r1_fastq%_R1*}_sorted_markdup.bam"

	_align "$align_cmd" "$r1_fastq" "$r2_fastq" "$reference" 2>"$log_file" \
		| samtools sort -@ "$(nproc)" -n - \
		| samtools fixmate -@ "$(nproc)" -m - - \
		| samtools sort -@ "$(nproc)" - \
		| samtools markdup -@ "$(nproc)" - "$bam_file"
}

run_on_all_fastq_files_in_dir() {
	local dir="$1"
	local align_cmd="hisat2"
	local reference="$HISAT2_MOUSE_INDEX"

	cd "$dir" || {
		echo "Error: $dir does not exist"
		exit 1
	}

	for fastq_file in ./*_R1.fastq.gz; do
		local r1_fastq="$fastq_file"
		local r2_fastq="${r1_fastq%_R1*}_R2.fastq.gz"
		_align2bam "$align_cmd" "$r1_fastq" "$r2_fastq" "$reference"
	done
}

run_on_all_fastq_files_in_dir "$1"
