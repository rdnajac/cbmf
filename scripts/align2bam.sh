#!/bin/bash
## Wrapper script for paired-end read alignment using HISAT2.
## Takes input and output directories as arguments.

set -euo pipefail

INPUT_FOLDER="20241029_Jurkat_PHF6-PHIP-KO_LQ/"
OUTPUT_FOLDER="$HOME/aligned"
mkdir -p "$OUTPUT_FOLDER"
THREADS="$(nproc)"
REFERENCE="/opt/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set"

cd "$INPUT_FOLDER" || {
	echo "Error: $INPUT_FOLDER does not exist"
	exit 1
}

for r1_fastq in ./*_R1*.fastq.gz; do
	r2_fastq="${r1_fastq/_R1/_R2}"

	if [[ ! -f "$r2_fastq" ]]; then
		echo "Warning: No matching R2 file found for $r1_fastq"
		continue
	fi

	# Extract sample name without _R1 and suffix
	sample_name="${r1_fastq##*/}"
	sample_name="${sample_name/_R1_*/}"
	bam_file="$OUTPUT_FOLDER/${sample_name}_sorted_markdup.bam"

	hisat2 -p "$THREADS" --dta -x "$REFERENCE" -1 "$r1_fastq" -2 "$r2_fastq" \
		| samtools sort -n -@ "$THREADS" - \
		| samtools fixmate -@ "$THREADS" -m - - \
		| samtools sort    -@ "$THREADS" - \
		| samtools markdup -@ "$THREADS" - "$bam_file"

	echo "Alignment completed for $r1_fastq and $r2_fastq. Output BAM: $bam_file"
done
