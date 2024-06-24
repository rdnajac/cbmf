#!/bin/bash

# Set variables for scripts
gtf="/home/ubuntu/ref/GRCm39/gencode.vM31.annotation.gtf"
projectDir="/home/ubuntu/export"
fastqDir="${projectDir}/00_fastq"
alignedDir="${projectDir}/aligned"
nCore=16

cd ${projectDir}
featureCounts \
	-a $gtf \
	-o DNMT3A_RHOA_PreLymphoma_APL_gencode.vM31.counts \
	-g 'gene_name' \
	-p -P -C -B \
	-T $nCore \
	${alignedDir}/*.bam
