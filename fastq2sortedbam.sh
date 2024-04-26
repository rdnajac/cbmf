#!/bin/bash
set -xeuo pipefail

THREADS=64
DIRECTORY="${1:-.}"
DIRECTORY="${DIRECTORY%/}"  # sanitize path
HUMANREF="/home/ubuntu/genomes/bowtie2/human/hg38"
MOUSEREF="/home/ubuntu/genomes/bowtie2/mouse/m39"
FASTQ_IDS=$(ls $DIRECTORY | grep _R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//')
IFS=$'\n' FASTQ_IDS=($FASTQ_IDS) # 

for sample in "${FASTQ_IDS[@]}"; do
    LOGFILE="${DIRECTORY}/${sample}_log.txt"
    SORTED_BAM="${DIRECTORY}/${sample}.sorted.bam"
    SAM_FILE="${DIRECTORY}/${sample}.sam"

    # Bowtie2 alignment with output written directly to the log file
    { time bowtie2 -p "$THREADS" -x "$MOUSEREF" \
                   -1 "${DIRECTORY}/${sample}_R1_001.fastq.gz" \
                   -2 "${DIRECTORY}/${sample}_R2_001.fastq.gz" \
                   -S "$SAM_FILE" 
    } 2>&1 | tee -a "$LOGFILE"

    # convert SAM to BAM, sort, and remove intermediate SAM file
    samtools view -buS "$SAM_FILE" | samtools sort -@ "$THREADS" -O bam -l 0 -T /tmp -o "$SORTED_BAM"
    rm -f "$SAM_FILE"
done
