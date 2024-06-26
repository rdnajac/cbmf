#!/bin/bash

# Declare associative arrays for genome and index paths
declare -A GENOME_RELEASES=(
	["human"]="GCA_000001405.15_GRCh38_full_analysis_set"
	["mouse"]="GCA_000001635.9_GRCm39_full_analysis_set"
)

declare -A GENOME_REFERENCES
declare -A HISAT2_INDEXES
declare -A BOWTIE_INDEXES

LOCAL_GENOMES_DIRECTORY="$HOME/genomes"

# Initialize paths
for organism in "${!GENOME_RELEASES[@]}"; do
	GENOME_REFERENCES[$organism]="${LOCAL_GENOMES_DIRECTORY}/${organism}/${GENOME_RELEASES[$organism]}"
	HISAT2_INDEXES[$organism]="${GENOME_REFERENCES[$organism]}.fna.hisat2_index"
	BOWTIE_INDEXES[$organism]="${GENOME_REFERENCES[$organism]}.fna.bowtie_index"
done

# Print all variables
_test() {
	echo "LOCAL_GENOMES_DIRECTORY: $LOCAL_GENOMES_DIRECTORY"
	for organism in "${!GENOME_RELEASES[@]}"; do
		echo "LATEST_MAJOR_${organism^^}_RELEASE: ${GENOME_RELEASES[$organism]}"
		echo "${organism^^}_REFERENCE: ${GENOME_REFERENCES[$organism]}"
		echo "HISAT2_${organism^^}_INDEX: ${HISAT2_INDEXES[$organism]}"
		echo "BOWTIE_${organism^^}_INDEX: ${BOWTIE_INDEXES[$organism]}"
	done
}

_export() {
	export LOCAL_GENOMES_DIRECTORY
	for organism in "${!GENOME_RELEASES[@]}"; do
		export "LATEST_MAJOR_${organism^^}_RELEASE=${GENOME_RELEASES[$organism]}"
		export "${organism^^}_REFERENCE=${GENOME_REFERENCES[$organism]}"
		export "HISAT2_${organism^^}_INDEX=${HISAT2_INDEXES[$organism]}"
		export "BOWTIE_${organism^^}_INDEX=${BOWTIE_INDEXES[$organism]}"
	done
}

_append_to_dot_bashrc() {
	{
		echo "export LOCAL_GENOMES_DIRECTORY=$LOCAL_GENOMES_DIRECTORY"
		for organism in "${!GENOME_RELEASES[@]}"; do
			echo "export LATEST_MAJOR_${organism^^}_RELEASE=${GENOME_RELEASES[$organism]}"
			echo "export ${organism^^}_REFERENCE=${GENOME_REFERENCES[$organism]}"
			echo "export HISAT2_${organism^^}_INDEX=${HISAT2_INDEXES[$organism]}"
			echo "export BOWTIE_${organism^^}_INDEX=${BOWTIE_INDEXES[$organism]}"
		done
	} >> "$HOME"/.bashrc
}

[[ "$1" = "test" ]] && _test && exit 0
[[ "$1" = "append" ]] && _append_to_dot_bashrc && exit 0
[[ "$1" = "export" ]] && _export && exit 0
[[ "$1" = "help" ]] && echo "Usage: $0 [test|append|export|help]" && exit 0
