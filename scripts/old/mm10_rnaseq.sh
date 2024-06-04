#!/bin/bash

for i in $(ls *_R1_001.fastq.gz)
do
f=${i%_R1_001.fastq.gz}
hisat2 -p 4 -x /Users/Robert/Desktop/index_chr19_and_others/hisat2_mm10/genome -1 ${f}_R1_001.fastq.gz -2 ${f}_R2_001.fastq.gz -S ${f}.sam
samtools view -b ${f}.sam > ${f}.test.bam
samtools sort -o ${f}.sorted.bam ${f}.test.bam
samtools index ${f}.sorted.bam
rm ${f}.sam 
rm ${f}.test.bam 

htseq-count -f bam ${f}.sorted.bam /Users/Robert/Desktop/index_chr19_and_others/mm10/gencode.vM24.annotation.gff3 -r pos > ${f}.counsHTSeq.csv

done

