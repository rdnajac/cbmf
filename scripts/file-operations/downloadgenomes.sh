#!/bin/bash

GENOMES_FOLDER=/opt/genomes
HUMAN_FOLDER=$GENOMES_FOLDER/human/GRCh38
MOUSE_FOLDER=$GENOMES_FOLDER/mouse/GRCm39

HUMAN_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
MOUSE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/"

# Create directories with sudo if they don't exist
sudo mkdir -p $HUMAN_FOLDER
sudo mkdir -p $MOUSE_FOLDER

# Download specific human genome files
echo "Downloading human genome files..."
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz" &
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz" &
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.fna.fai" &
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz" &
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.fna.hisat2_index.tar.gz" &
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz" &
sudo wget -c -P $HUMAN_FOLDER "$HUMAN_URL/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz" &

# Download entire mouse genome folder
echo "Downloading entire mouse genome folder..."
sudo wget -r -np -nH --cut-dirs=7 -P $MOUSE_FOLDER "$MOUSE_URL"

echo "Download complete."
