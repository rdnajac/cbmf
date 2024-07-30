#!/bin/bash
#
# This file acts as a wrapper around the StringTie program.
# Run these commands after aligning the RNA-seq reads to the genome.
#
# Input: sorted bam files
# Output: GTF files
set -euo pipefail
IFS=$'\n\t'

set -x

# mouse
REF_GFF="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff"
REF_GTF="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gtf"

MERGED_GTF="stringtie_merged.gtf"

# Get the number of processors and save it to a variable
NPROC=$(nproc)

cleanup() {
	mv -v ./gffcmp* ./gffcompare_output
	mv -v ./*.gtf ./stringtie_output
	rm -rfv tmp*
}

prepare_for_export() {
	local prefix="$(basename "$PWD")"

	zip -rv "$prefix"_stringtie.zip stringtie_output gffcompare_output ballgown_input
	rm -frv stringtie_output gffcompare_output ballgown_input
}

# main function
main() {
	local directory="${1:-.}"
	local species="${2:-mouse}"

	cd "$directory" || {
		echo "Directory not found."
		exit 1
	}

	if [[ "$species" != "mouse" ]]; then
		echo "Species not supported."
		exit 2
	fi

	# output folders
	mkdir -vp ballgown_input
	mkdir -vp stringtie_output
	mkdir -vp gffcompare_output

	reference_dir="$HOME/genomes/mouse"
	reference_gff="$reference_dir/$REF_GFF"
	reference_gtf="$reference_dir/$REF_GTF"

	# Run StringTie on each sample to assemble and quantify expressed genes and transcripts
	for bamfile in ./*.bam; do
		stringtie -p "$NPROC" -G "$reference_gff" "$bamfile" > "${bamfile%.bam}_ref.gtf" &
	done
	echo "Waiting for transcripts to be assembled and quantified..."
	wait
	echo "Transcripts assembled and quantified."

	# 2. merge GTF files
	echo "Merging GTF files..."
	stringtie --merge -p "$NPROC" -G "$reference_gff" ./*_ref.gtf > "$MERGED_GTF"
	echo "Merged GTF file created: $MERGED_GTF"

	# gffcompare to compare the merged GTF file to the reference annotation
	# this can run in parallel
	{
		echo "Running gffcompare..."
		gffcompare "$reference_gff" "$MERGED_GTF" -V 2> gffcmp.log
	} &

	# 4. generate Ballgown input files
	echo "Generating Ballgown input files..."
	for bamfile in ./*.bam; do
		local id="${bamfile%.bam}"
		local folder="ballgown_input/$id"
		mkdir -vp "$folder"
		stringtie -p "$NPROC" -G "$MERGED_GTF" -e -b "$folder" -o "${id}_stringtie.gtf" "$bamfile" &
	done
	wait
	echo "All samples processed with StringTie using merged GTF."
	cleanup
	prepare_for_export
}

main "$@"
