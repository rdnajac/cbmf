#!/bin/bash
#
# This script is a pipeline to align paired-end reads to a reference (mouse) genome

# Enforce "strict" error handling
set -euo pipefail

set -x

export GENOMES="/home/ubuntu/genomes"
export MOUSEREF="${GENOMES}/mouse/GCA_000001635.9_GRCm39_full_analysis_set"
export HISAT2_MOUSE_INDEX="${MOUSEREF}.fna.hisat2_index"
export BOWTIE_MOUSE_INDEX="${MOUSEREF}.fna.bowtie_index"
export BOWTIE_HUMAN_INDEX="${GENOMES}/human/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"

export HUMAN_BOWTIE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz"

# These flags are used by both aligners
COMMON_FLAGS="-p $(nproc) --mm"

# BOWTIE2_FLAGS="${COMMON_FLAGS} -x ${BOWTIE_MOUSE_INDEX}"
BOWTIE2_FLAGS="${COMMON_FLAGS} -x ${BOWTIE_HUMAN_INDEX}"

# use flag -dta for --downstream-transcriptome-assembly
# HISAT2_FLAGS="${COMMON_FLAGS} --dta -x ${HISAT2_MOUSE_INDEX}"

fastq2bam()
{
	local id="$1"
	local suffix="${2:-_001.fastq.gz}"
	local r1="${id}_R1$suffix"
	local r2="${id}_R2$suffix"
	local log_file="${id}_align.log"
	local bam_file="${id}_sorted_markdup.bam"

	# hisat2 "$HISAT2_FLAGS" -1 "$r1" -2 "$r2" --met-stderr 2> >(tee "$log_file" >&2) \
	bowtie2 "$BOWTIE2_FLAGS" -1 "$r1" -2 "$r2" 2> >(tee "$log_file" >&2) \
		| samtools sort -@ "$(nproc)" -n - \
		| samtools fixmate -@ "$(nproc)" -m - - \
		| samtools sort -@ "$(nproc)" - \
		| samtools markdup -@ "$(nproc)" - "$bam_file"

	# optionally, write protect the bam file
	chmod a-w "$bam_file"
}

cd "$1" || exit 1

for f in ./*R1_001.fastq.gz; do
	id="${f%%_R1_001.fastq.gz}"
	fastq2bam "$id"
done
