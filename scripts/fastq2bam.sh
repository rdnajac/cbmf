#!/bin/bash
# run this script in the folder to process
set -euox pipefail
THREADS=$(nproc)
SUFFIX="_001.fastq.gz"
REFERENCE_GENOME="$HUMANREF" # change this
FASTA="/home/ubuntu/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"


align2sortedbam() {
  # ~sage: align2sortedbam <ref> <id>
  local ref=$1
  local id=$2
  local r1_file="${id}_R1${SUFFIX}"
  local r2_file="${id}_R2${SUFFIX}"
  local logfile="${id}.log"
  local bamfile="${id}.bam"

    # time bowtie2 -p "${THREADS}" --mm -x "${REFERENCE_GENOME}" \
    #   -1 "${r1_file}" -2 "${r2_file}" 2>> "${id}.log" | \
    # time samtools sort -@ "$THREADS" -o "${id}.bam" 2>> "${id}.log"
    time {
    bowtie2 -p "${THREADS}" --mm -x "${ref}" -1 "${r1_file}" -2 "${r2_file}" 2>> "${logfile}" | \
      samtools sort -@ "$THREADS" -o "${bamfile}" 2>> "${logfile}"
    }
}

#SAMPLES=$(~/index.sh .)
#make sure its an array we might get space or newline separated
declare -a SAMPLES=($(~/index.sh .))

for sample in "${SAMPLES[@]}"; do
  # change this to the reference genome you want to use
  # TODO: this is lazy, fix it
  align2sortedbam $HUMANREF $sample
done
