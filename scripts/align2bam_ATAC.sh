#!/bin/bash

set -euo pipefail

INPUT_FOLDER="$HOME/ATACseq"
OUTPUT_FOLDER="$HOME/aligned"
mkdir -p "$OUTPUT_FOLDER"
THREADS="$(nproc)"
moREFERENCE="/opt/genomes/mouse/GRCm39/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index"
huREFERENCE="/opt/genomes/human/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"

human_samples=("MyLa" "MJ" "HuT-78" "T8ML-1" "OCI-Ly12" "PDX-78024")
mouse_samples=("DNMT3A" "FYN-TRAF3IP2" "Tet2Rhoa" "Vav1-Myo1f" "muCD4" "muTFH")

cd "$INPUT_FOLDER" || {
	echo "Error: $INPUT_FOLDER does not exist"
	exit 1
}

fastq2bam() { 
	sampleID=$1
	refGenome=$2
	r1_fastq="${sampleID}_R1_001.fastq.gz"
	r2_fastq="${sampleID}_R2_001.fastq.gz"
	bam_file="${OUTPUT_FOLDER}/${sampleID}_sorted_markdup.bam"

	bowtie2 -p "$THREADS" -x "$refGenome" -1 "$r1_fastq" -2 "$r2_fastq" \
		| samtools sort -n -@ "$THREADS" - \
		| samtools fixmate -@ "$THREADS" -m - - \
		| samtools sort    -@ "$THREADS" - \
		| samtools markdup -@ "$THREADS" - "$bam_file"

	echo "Alignment completed for ${sampleID}. Output BAM: $bam_file"
}

for mousesample in "${mouse_samples[@]}"; do
	fastq2bam "$mousesample" "$moREFERENCE"
done

for humansample in "${human_samples[@]}"; do
	fastq2bam "$humansample" "$huREFERENCE"
done
