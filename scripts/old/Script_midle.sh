
#!/bin/bash

bowtie2 --threads 5 --very-sensitive --no-unal --no-mixed  --phred33 -X 700 -x /Users/Robert/Desktop/mm10/GRCm38.6.plusribosomic -U RA678-H3K4me1_S1_L001_R1_001.fastq.gz,RA678-H3K4me1_S1_L002_R1_001.fastq.gz,RA678-H3K4me1_S1_L003_R1_001.fastq.gz,RA678-H3K4me1_S1_L004_R1_001.fastq.gz  -S H3K4me1.sam
samtools view -@ 12 -b H3K4me1.sam > H3K4me1.bam
samtools sort -@ 12  H3K4me1.bam -o H3K4me1.s.bam
samtools index H3K4me1.s.bam
bamCoverage --bam H3K4me1.s.bam -o H3K4me1.bw -p 16 --binSize 20 --normalizeUsing RPKM --smoothLength 30 --effectiveGenomeSize 2864785220

bowtie2 --threads 5 --very-sensitive --no-unal --no-mixed  --phred33 -X 700 -x /Users/Robert/Desktop/mm10/GRCm38.6.plusribosomic -U RA678-H3K27ac_S2_L001_R1_001.fastq.gz,RA678-H3K27ac_S2_L002_R1_001.fastq.gz,RA678-H3K27ac_S2_L003_R1_001.fastq.gz,RA678-H3K27ac_S2_L004_R1_001.fastq.gz -S H3K27ac.sam
samtools view -@ 12 -b H3K27ac.sam > H3K27ac.bam
samtools sort -@ 12  H3K27ac.bam -o H3K27ac.s.bam
samtools index H3K27ac.s.bam
bamCoverage --bam H3K27ac.s.bam -o H3K27ac.bw -p 16 --binSize 20 --normalizeUsing RPKM --smoothLength 30 --effectiveGenomeSize 2864785220
