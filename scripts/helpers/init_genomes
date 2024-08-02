#!/bin/bash
#
## This script installs the reference genomes for the pipelines
## and sets up the necessary directories and environment variables.

set -euo pipefail
set -x

# create the genomes directory
GENOMES_DIRECTORY="$HOME/genomes"
mkdir -vp "$GENOMES_DIRECTORY"
cd "$GENOMES_DIRECTORY" || {
	echo "Error: $GENOMES_DIRECTORY does not exist"
	exit 1
}

# download the reference genomes
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39

NCBI="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
MOUSE="000/001/635/GCA_000001635.9_GRCm39/"
HUMAN="000/001/405/GCA_000001405.15_GRCh38/"

echo "full mouse url: $NCBI/$MOUSE"
echo "full human url: $NCBI/$HUMAN"

# download the reference entire folder
# use the -e robots=off flag to ignore the robots.txt file
# which normally prevents web crawlers from downloading the entire site 
#
# use the --cut-dirs flag to flatten the directory structure so that
# the files are downloaded to the current directory
# wget --recursive -e robots=off --no-parent --cut-dirs=9 "$NCBI/$MOUSE/"

download_reference_genome()
{
	local url="$1"
	wget --recursive -e robots=off --no-parent --cut-dirs=9 "$url"
}

# download_reference_genome "$NCBI/$MOUSE/seqs_for_alignment_pipelines.ucsc_ids/"
# download_reference_genome "$NCBI/$HUMAN/seqs_for_alignment_pipelines.ucsc_ids/"

# notably, we get pre-built indexes for bwa, HISAT2, Bowtie2
# as well as annotation files in .gff and .gtf formats

