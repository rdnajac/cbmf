#!/bin/bash

set -euo pipefail   # enforce strict error handling
set -x

THREADS=8
REF="/home/ubuntu/rnaseq/genomes/mouse2/GCA_000001635.9_GRCm39_full_analysis_set.fna.hisat2_index"
ANNOTATION="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gtf"

SAMPLES=(
  DMSO{1..3}
  Fingolimod{1..3}
  Ozanimod{1..3}
  Ponesimod{1..3}
)

fastq2bam() {
  local id="$1"
  local r1="${id}_R1_001.fastq.gz"
  local r2="${id}_R2_001.fastq.gz"
  local log_file="${id}.log"
  local bam_file="${id}_sorted_markdup.bam"

  # Align the reads and process with samtools
  hisat2 --time --threads "$THREADS" --mm -x "$REF" -1 "$r1" -2 "$r2" --met-stderr 2> "$log_file" \
  | samtools sort -n -@ "$THREADS" - \
  | samtools fixmate -m -@ "$THREADS" - - \
  | samtools sort -@ "$THREADS" - \
  | samtools markdup -@ "$THREADS" - "$bam_file"
}

# for id in "${SAMPLES[@]}"; do
#   fastq2bam "$id" &
# done

STRINGTIE="stringtie -p 8"
for id in "${SAMPLES[@]}"; do
  ${STRINGTIE} -G "${ANNOTATION}" -o "${id}.gtf" "${id}_sorted_markdup.bam" &
done

# Merge all the GTF files
#${STRINGTIE} --merge -G "${ANNOTATION}" -o stringtie_merged.gtf *.gtf

# for id in "${SAMPLES[@]}"; do
#   # stringtie -e -B -p "$THREADS" -G stringtie_merged.gtf -o "${id}.gtf" "${id}_sorted_markdup.bam" &
#   # properly name hte outpt
#   stringtie -e -B -p "$THREADS" -G stringtie_merged.gtf -o "${id}_stringtie.gtf" "${id}_sorted_markdup.bam" &
# done
