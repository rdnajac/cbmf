#!/bin/bash
#
# This file acts as a wrapper around the StringTie program.
# Run these commands after aligning the RNA-seq reads to the genome.
#
# Input: sorted bam files
# Output: GTF files
set -euo pipefail
set -x

GENOMES="/home/ubuntu/genomes"
MOUSEREF="${GENOMES}/mouse/GCA_000001635.9_GRCm39_full_analysis_set"
MREF_GTF="${MOUSEREF}.refseq_annotation.gtf"
MREF_GFF="${MOUSEREF}.refseq_annotation.gff"

## Run StringTie on a single sample to assemble and quantify expressed genes and transcripts
run_stringtie()
{
	local reference="$1"
	local input_bam="$2"
	local output_gtf="$3"

	stringtie -p "$(nproc)" -G "$reference" "$input_bam" >"$output_gtf"
}

## Run StringTie on each sample in the current directory
do_run_stringtie()
{
	local reference="$1"
	local input_folder="${2:-.}"
	local output_folder="${3:-.}"

	for f in "$input_folder"/*.bam; do
		run_stringtie "$reference" "$f" "${output_folder}/$(basename "${f%.bam}").gtf" &
	done
	wait
	echo "All samples processed with StringTie."
}

# File to store the merged GTF information
MERGED_GTF="stringtie_merged.gtf"

# Merge GTF files
# @sideeffect creates a merged GTF file
merge_gtfs()
{
	local input_folder="${1:-.}"

	echo "Merging GTF files using GFF"
	stringtie --merge -p "$(nproc)" -G "$MREF_GFF" "$input_folder"/*.gtf >"$MERGED_GTF"
	echo "GTF files merged into $MERGED_GTF."
}

# Compare GFF files
do_gffcompare()
{
	local reference="$1"
	local merged_gtf="$2"
	gffcompare -r "$reference" -G -o merged "$merged_gtf"
}

# Generate Ballgown input files
# @brief Process samples with merged GTF
generate_ballgown_inputs()
{
	local reference="$1"
	local input_folder="${2:-.}"

	for f in "$input_folder"/*.bam; do
		local id="${f%.bam}"
		local folder="ballgown_input/$id"
		mkdir -vp "$folder"
		echo "stringtie -p $(nproc) -G $reference -e -b $folder -o ${id}_stringtie.gtf $f"
		stringtie -p "$(nproc)" -G "$reference" -e -b "$folder" -o "${id}_stringtie.gtf" "$f" &
	done
	wait
	echo "All samples processed with StringTie using merged GTF."
}

# 1. run StringTie on each sample in pwd

FILES=(
	./Tet2RhoaG17V_24hAGX51_1_sorted_markdup.bam
	./Tet2RhoaG17V_24hAGX51_2_sorted_markdup.bam
	./Tet2RhoaG17V_24hAGX51_3_sorted_markdup.bam
	./Tet2RhoaG17V_24hDMSO_1_sorted_markdup.bam
	./Tet2RhoaG17V_24hDMSO_2_sorted_markdup.bam
	./Tet2RhoaG17V_24hDMSO_3_sorted_markdup.bam
)

#same but uses our FILES
do_run_stringtie_on_files()
{
	local reference="$1"

	for file in "${FILES[@]}"; do
		run_stringtie "$reference" "$file" "${file%.bam}.gtf" &
		echo "StringTie on $file."
	done
	wait
	echo "All samples processed with StringTie."
}

# do_run_stringtie "$MREF_GFF"

# 2. merge GTF files
# merge_gtfs ~/FASTQ/gtf

# 3. compare GFF files
# cd ~/FASTQ/gtf
# do_gffcompare "$MREF_GFF" "$MERGED_GTF"

# 4. generate Ballgown input files
cd ~/FASTQ/bam
generate_ballgown_inputs "$MREF_GFF"

## notes
# the default output foler contains both the .gtf and _stringtie.gtf files
