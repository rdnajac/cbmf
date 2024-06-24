#!/bin/bash

# Set variables for scripts
gtf="/mnt/data/shares/ref/STAR/GRCh38p12/gencode.v31.annotation.gtf"
projectDir="/mnt/data/bobby/scratch/RNAseq/20220517_MYLA_RomiAfaDimSyn"
fastqDir="${projectDir}/00_fastq"
alignedDir="${projectDir}/aligned"
nCore=12

cd ${projectDir}
featureCounts \
	-a $gtf \
	-o MYLA_RomiAfaDimSyn_GRCh38p12_gencodev31.counts \
	-p -P -C -B \
	-T $nCore \
	${alignedDir}/*.bam
