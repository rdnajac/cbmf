#!/bin/bash

set -euo pipefail   # enforce strict error handling
set -x

THREADS="$(nproc)"

HISAT_REF_INDEX="/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.fna.hisat2_index"
REFSEQ_ANNOTATION_GTF="/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gtf"
REFSEQ_ANNOTATION_GFF="/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff"

# Align the reads and process with samtools
# The first line is variable depending on which software used for alignment
# the samtools pipeline should be the same for all aligners that output SAM
fastq2bam() {
  local id="$1"
  local r1="${id}_R1_001.fastq.gz"
  local r2="${id}_R2_001.fastq.gz"
  local log_file="${id}_hisat2.log"
  # local bam_file="${id}_sorted_markdup.bam"
  local bam_file="${id}.bam"

  hisat2 --mm --threads "$THREADS" -x "$HISAT_REF_INDEX" -1 "$r1" -2 "$r2" --met-stderr 2> >(tee "$log_file" >&2) \
  | samtools sort    -@ "$THREADS" -n -   \
  | samtools fixmate -@ "$THREADS" -m - - \
  | samtools sort    -@ "$THREADS"    -   \
  | samtools markdup -@ "$THREADS"    - "$bam_file"
}

for r1 in *_R1_001.fastq.gz; do
  id="${r1%%_R1_001.fastq.gz}"
  fastq2bam "$id" &
  # wait for 3 jobs to finish
  if [[ $(jobs -p | wc -l) -ge 3 ]]; then
    wait
  fi
done

# STRINGTIE="stringtie -p 8"
# for id in "${SAMPLES[@]}"; do
#   ${STRINGTIE} -G "${ANNOTATION}" -o "${id}.gtf" "${id}_sorted_markdup.bam" &
# done
#
# one-liner to run stringtie on all samples with nproc threads (use ./*.bam to avoid issues with the order of the files)
for f in ./*.bam; do stringtie -p "$THREADS" -G "$REFSEQ_ANNOTATION_GTF" -o "${f%.bam}.gtf" "$f" & done
# use gtf because it is more 'gene-centric' than gff
for f in ./*.bam; do stringtie -p 8 -G "/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff" -l ${f} -o "${f%.bam}.gtf" "$f" & done

# merge those files
echo $(ls *.gtf) > merge_list.txt && cat merge_list.txt
# do it again but use the newlines from ls
ls *.gtf > merge_list.txt && cat merge_list.txt
stringtie --merge -G "$REFSEQ_ANNOTATION_GTF" -o stringtie_merged.gtf ./*.gtf
stringtie --merge -p 32 -G /home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff -o stringtie_merged.gtf merge_list.txt

# now re-run stringtie with the merged file on each sample
for f in ./*.bam; do stringtie -eB -p 32 -G stringtie_merged.gtf -o "${f%.bam}_stringtie.gtf" "$f" & done

# for id in "${SAMPLES[@]}"; do

# Merge all the GTF files
#${STRINGTIE} --merge -G "${ANNOTATION}" -o stringtie_merged.gtf *.gtf

# for id in "${SAMPLES[@]}"; do
#   # stringtie -e -B -p "$THREADS" -G stringtie_merged.gtf -o "${id}.gtf" "${id}_sorted_markdup.bam" &
#   # properly name hte outpt
#   stringtie -e -B -p "$THREADS" -G stringtie_merged.gtf -o "${id}_stringtie.gtf" "${id}_sorted_markdup.bam" &
# done
#
# write protect all *.fastq.gz and *.bam files
chmod a-w ./*.fastq.gz ./*.bam
