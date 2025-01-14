#!/bin/bash

for i in $(ls *_L001_R1_001.fastq.gz)
do
f=${i%_L001_R1_001.fastq.gz}
bowtie2 --no-unal --threads 8 --local --very-sensitive-local  --no-unal --no-mixed --no-discordant  --phred33 -X 700 -t -x /Users/Robert/Desktop/mm10/GRCm38.6.plusribosomic -U ${f}_L001_R1_001.fastq.gz,${f}_L002_R1_001.fastq.gz,${f}_L003_R1_001.fastq.gz,${f}_L004_R1_001.fastq.gz -S ${f}.sam
samtools view -@ 11 -b ${f}.sam> ${f}.bam
samtools sort -@ 11 -m 8G -o ${f}.sorted.bam ${f}.bam
samtools index ${f}.sorted.bam
bamCoverage -b ${f}.sorted.bam -o ${f}.sorted.normRPKM.bw --normalizeUsing RPKM
done