#!/bin/bash
#
## Convert BAM to CRAM format

set -euo pipefail

# url for human download:
# url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

# download human genome
# cd /home/ubuntu/genomes/human && curl -O -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

# gunzip it and bgzip it
# gunzip -c GCA_000001405.15_GRCh38_full_analysis_set.fna.gz | bgzip -c > GCA_000001405.15_GRCh38_full_analysis_set.fna.bgz
# index the genome

bam2cram()
{
	local bam_file="$1"
	local id="$(basename "$bam_file" .bam)"
	local cram_file="${id}.cram"
	local mouse_fai="/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.fna.bgz"
	local human_fai="/home/ubuntu/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set.fna.bgz"

	samtools view -@"$(nproc)" --cram -T "$human_fai" "$bam_file" > "$cram_file"

	# # optionally, write protect the cram file
	# chmod a-w "$cram_file"

	# # remove the bam file
	# chmod a+w "$bam_file" && rm "$bam_file"
}

cram2bam()
{
	local cram_file="$1"
	local id="$(basename "$cram_file" .cram)"
	local bam_file="${id}.bam"

	# $if bamfile aready exists, move it to a backup file
	if [ -f "$bam_file" ]; then
		mv "$bam_file" "${bam_file}.bak"
	fi

	samtools view -@"$(nproc)" --bam "$cram_file" > "$bam_file"

	# # optionally, write protect the bam file
	# chmod a-w "$bam_file"

	# # remove the cram file
	# chmod a+w "$cram_file" && rm "$cram_file"
}

# switch case for suffix o function call
if [ $# -eq 0 ]; then
	echo "Usage: $0 <bam_file>"
	exit 1
fi

case "$1" in
*.bam) bam2cram "$1" ;;
*.cram) cram2bam "$1" ;;
*) echo "Unknown file type: $1" ;;
esac
