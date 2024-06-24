#!/bin/bash
#
# This script is a pipeline to align paired-end reads to a reference (mouse) genome

# Enforce "strict" error handling
set -euo pipefail

export GENOMES="/home/ubuntu/genomes"
export MOUSEREF="${GENOMES}/mouse/GCA_000001635.9_GRCm39_full_analysis_set"
export HISAT_MOUSE_INDEX="${MOUSEREF}.fna.hisat2_index"
export MREF_GTF="${MOUSEREF}.refseq_annotation.gtf"
export MREF_GFF="${MOUSEREF}.refseq_annotation.gff"

## Parse command line options to set flags with string building

# These flags are used by both aligners
COMMON_FLAGS="-p $(nproc) --mm"
# BOWTIE2_FLAGS="${COMMON_FLAGS} --local"
HISAT2_FLAGS="${COMMON_FLAGS} --dta -x ${HISAT_MOUSE_INDEX}"

fastq2bam()
{
	local id="$1"
	local suffix="$2" # .fastq.gz or _001.fastq.gz
	local r1="${id}_R1$suffix"
	local r2="${id}_R2$suffix"
	local log_file="${id}_hisat2.log"
	local bam_file="${id}_sorted_markdup.bam"

	hisat2 "$HISAT2_FLAGS" -1 "$r1" -2 "$r2" --met-stderr 2> >(tee "$log_file" >&2) \
		| samtools sort -@ "$(nproc)" -n - \
		| samtools fixmate -@ "$(nproc)" -m - - \
		| samtools sort -@ "$(nproc)" - \
		| samtools markdup -@ "$(nproc)" - "$bam_file"

	# Optionally, remove duplicates

	# write protect the bam file
	chmod a-w "$bam_file"
}

cd FASTQ || exit

# make a logfile for this run

# redirect stdout and stderr to the logfile
for f in *R1.fastq.gz; do
	id="${f%%_R1.fastq.gz}"
	fastq2bam "$id" ".fastq.gz"
done
