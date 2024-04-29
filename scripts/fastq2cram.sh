#!/bin/bash

set -euox pipefail
THREADS=$(nproc)

readonly SUFFIX="_001.fastq.gz"

# Function to align a pair of reads to a reference genome and convert output directly to BAM format
fastq2bam() {
    local dir=$1
    local id=$2
    local ref=$3

    local r1_file="${id}_R1${SUFFIX}"
    local r2_file="${id}_R2${SUFFIX}"
    local outfile="${dir}/${id}.bam"
    local logfile="${dir}/${id}.log"

    bowtie2 -p "$THREADS" --mm -x "$ref" -1 "$r1_file" -2 "$r2_file" 2>> "$logfile" | \
    samtools sort -@ "$THREADS" -o "$outfile" 2>> "$logfile"
}

export hmFASTA="/home/ubuntu/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"

bam2cram() {
    local bam_file="$1"
    local cram_file="${bam_file%.bam}.cram"
    local log_file="${cram_file%.cram}.log"

    time samtools view -T "${hmFASTA}" -C -@ "$THREADS" -o "$cram_file" "$bam_file" 2>> "$log_file"

    echo "Conversion completed: $bam_file to $cram_file"
}

cram2bam() {
    local cram_file="$1"
    local bam_file="${cram_file%.cram}.bam"
    local log_file="${bam_file%.bam}.log"

    time samtools view -T "${hmFASTA}" -b -@ "$THREADS" -o "$bam_file" "$cram_file" 2>> "$log_file"

    echo "Conversion completed: $cram_file to $bam_file"
}

samtools view -T "${hmFASTA}"
