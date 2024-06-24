#!/bin/bash
#
## Convert BAM to CRAM format

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
	local mouse_fai="/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.fna.fai"

	# samtools view -@"$(nproc)" --cram -T /home/ubuntu/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set.fna.bgz "$bam_file" >"$cram_file"
	samtools view -@"$(nproc)" --cram -T "$mouse_fai" "$bam_file" >"$cram_file"

	# optionally, write protect the cram file
	chmod a-w "$cram_file"

	# remove the bam file
	chmod a+w "$bam_file" && rm "$bam_file"
}

cd "$1" || exit 1
for f in ./*.bam; do
	bam2cram "$f"
done
