#!/bin/bash -
# FILE: bam2cram.sh
# USAGE: ./bam2cram.sh
# DESCRIPTION:
# AUTHOR: Ryan D. Najac, ryan.najac@columbia.com
# ORGANIZATION: Columbia University
# CREATED: 2024-05-09 1715284585
# REVISION:  ---

set -euxo pipefail

readonly FNA="/home/ubuntu/download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz"
readonly MAX_THREADS="$(nproc)"

convert() {
  local file="$1"
  case "${file##*.}" in
    bam)
      samtools view -@ "$MAX_THREADS" -C -T "$FNA" -o "${file%.bam}.cram" "$file"
      echo "Converted $file to CRAM format"
      ;;
    cram)
      samtools view -@ "$MAX_THREADS" -b -o "${file%.cram}.bam" "$file"
      echo "Converted $file to BAM format"
      ;;
    *)
      echo "Error: Unsupported file type."
      exit 1
      ;;
  esac
}

for file in "$@"; do
  convert "$file"
done

