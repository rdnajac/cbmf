#!/bin/bash
#
## This script runs featureCounts on the BAM files in the input directory

# mouse
# REF_GFF="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff"
export REF_GTF="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gtf"

# other variavles
NUM_THREADS=$(nproc)

main() {
	local dir="$1"
	local fname="$2"


	cd "$dir" || {
		echo "Directory not found."
		exit 1
	}

	featureCounts -T "$NUM_THREADS" --verbose \
	        -t exon -g gene_id --countReadPairs \
		-a "${HOME}/genomes/mouse/${REF_GTF}" \
		-p -P -C -B -o "${fname}.tsv" ./*.bam
}


featureCounts -T $(nproc) --verbose -t exon -g gene_id --countReadPairs -a ${HOME}/genomes/mouse/${REF_GTF} -p -P -C -B -o fin_dmso.tsv ./F*.bam ./DM*.bam
