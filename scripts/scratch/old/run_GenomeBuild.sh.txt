#!/bin/bash

# 2021-05-28
# Script to generate genome

samples_file="/home/ubuntu/scratch/scripts/samples.txt"
genome_dir="/home/ubuntu/scratch/genome"

n_cores=`grep -c ^processor /proc/cpuinfo`

STAR \
	--runThreadN ${n_cores} \
	--runMode genomeGenerate \
	--genomeDir ${genome_dir}/GRCh38_150 \
	--genomeFastaFiles ${genome_dir}/GRCh38.primary_assembly.genome.fa \
	--sjdbGTFfile ${genome_dir}/gencode.v38.annotation.gtf \
	--sjdbOverhang 149
