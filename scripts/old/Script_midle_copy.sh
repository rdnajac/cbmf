
#!/bin/bash

bowtie2 --threads 8 --very-sensitive --no-unal --no-mixed  --phred33 -X 700 -x /Users/Robert/Desktop/mm10/GRCm38.6.plusribosomic -U Undetermined_S0_L001_R1_001.fastq.gz,Undetermined_S0_L002_R1_001.fastq.gz,Undetermined_S0_L003_R1_001.fastq.gz,Undetermined_S0_L004_R1_001.fastq.gz  -S unindex.sam
samtools view -@ 12 -b unindex.sam > unindex.bam
samtools sort -@ 12  unindex.bam -o unindex.s.bam
samtools index unindex.s.bam
bamCoverage --bam H3K4me1.s.bam -o H3K4me1.bw -p 16 --normalizeUsing RPKM 
bamCoverage --bam unindex.s.bam -o unindex.bw -p 16 --normalizeUsing RPKM 
bamCoverage --bam H3K27ac.s.bam -o H3K27ac.bw -p 16 --normalizeUsing RPKM 