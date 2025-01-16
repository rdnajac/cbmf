#!/bin/bash

GENOMES_FOLDER=/opt/genomes
MOUSE_FOLDER=$GENOMES_FOLDER/mouse/GRCm39

MOUSE_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/

# Create mouse folder if it doesn't exist
sudo mkdir -p $MOUSE_FOLDER

# Download mouse genome files in parallel using MOUSE_URL
echo "Downloading mouse genome files..."

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index.tar.gz" &

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.fna.bwa_index.tar.gz" &

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.fna.fai" &

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.fna.gz" &

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.fna.hisat2_index.tar.gz" &

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff.gz" &

sudo wget --no-check-certificate --directory-prefix=$MOUSE_FOLDER \
  "${MOUSE_URL}GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gtf.gz" &

# Wait for all background processes to finish
wait

echo "Download complete."
