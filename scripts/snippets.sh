q
#!/bin/bash

#================================================================
## Description: Utility functions and example scripts for RNAseq analysis
## Author: Ryan D. Najac
## Last modified: 2024-05-29
#================================================================

#======================= GLOBAL SETTINGS ================================
set -euxo pipefail   # enforce strict error handling

#======================= INSTALLATION FUNCTIONS =========================

# Function to install software from GitHub and add to PATH
install_and_add_to_path() {
  (
    git clone "$1" && cd "$(basename "$1" .git)" && make -j "$(nproc)"
    export PATH="$PATH:$(pwd)"
  )
}

# Install HISAT2
install_hisat2() {
  install_and_add_to_path https://github.com/DaehwanKimLab/hisat2.git
} q

# Install StringTie
install_stringtie() {
  install_and_add_to_path https://github.com/gpertea/stringtie.git
}

# Install STAR
install_star() {
  install_and_add_to_path https://github.com/alexdobin/STAR.git
}

# Install Ballgown in R
install_ballgown() {
  Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")'
  Rscript -e 'BiocManager::install("ballgown")'
}

#======================= CONVERSION UTILITIES ===========================

# Aliases for converting sample read files
# `source snippets.sh` to add them to your environment

alias fastq_to_fasta="sed 'N;x;N;N;x;s/@>/'"
alias paired_to_tab5="paste <(sed 'N;x;N;g;N;s/\n/\t/g' reads_1.fq) <(sed -n 'n;h;n;g;N;s/\n/\t/g;p' reads_2.fq) > reads_12.tab5"
alias paired_to_tab6="paste <(sed 'N;x;N;g;N;s/\n/\t/g' reads_1.fq) <(sed 'N;x;N;g;N;s/\n/\t/g' reads_2.fq) > reads_12.tab6"
alias paired_to_interleaved="paste -d'\n' <(sed 'N;N;N;s/\n/\t/g' reads_1.fq) <(sed 'N;N;N;s/\n/\t/g' reads_2.fq) | tr '\t' '\n' > reads_12.fq"

#======================= BAM TO CRAM CONVERSION ==========================

#!/bin/bash -
# FILE: bam2cram.sh
# USAGE: ./bam2cram.sh
# DESCRIPTION: Convert BAM files to CRAM format and vice versa
# AUTHOR: Ryan D. Najac
# CREATED: 2024-05-09

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

#======================= FASTQ TO BAM PIPELINE ==========================

#!/bin/bash
# Description: Fastq to aligned BAM pipeline
# Author: Ryan D. Najac
# Last modified: 2024-05-03

readonly THREADS=$(nproc)
readonly MOUSEREF="download/GCA_000001635.9_GRCm39_genomic.fna.gz"
readonly HUMANREF="/home/ubuntu/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"

get_fastq_pairs() {
  local dir=${1:-$(pwd)}
  local suffix="_001.fastq.gz"
  declare -a files=()

  while IFS= read -r r1; do
    sample_id=$(basename "$r1" "_R1$suffix")
    local r2="$dir/${sample_id}_R2$suffix"
    [[ -f "$r2" ]] && files+=("${dir}/${sample_id}")
  done < <(find "$dir" -maxdepth 1 -type f -regex ".*_R[12]$suffix")

  printf "%s\n" "${files[@]}"
}

fastq2bam() {
  local id="$1"
  local r1="${id}_R1_001.fastq.gz"
  local r2="${id}_R2_001.fastq.gz"
  local log_file="${id}.log"

  log() { echo "$(date +%H%M%S) $*" >> "$log_file"; }

  log "start $id"
  {
    bowtie2 --time --threads "$THREADS" --mm -x "$MOUSEREF" -1 "$r1" -2 "$r2"       | samtools sort -n -@ "$THREADS" -       | samtools fixmate -m -@ "$THREADS" - -       | samtools sort -@ "$THREADS" -       | samtools markdup -@ "$THREADS" - "${id}.bam"
  } 2> "$log_file"
  log "finish $id"
}

[[ -d "$1" ]] || { echo "Error: Directory '$1' not found."; exit 1; }

declare -a fastq_files=()
fastq_files=($(get_fastq_pairs "$1"))

for id in "${fastq_files[@]}"; do
  fastq2bam "$id"
done

#======================= FASTQC MULTIPLE FILES ==========================

#!/bin/bash
# Description: Script to run FastQC on various types of sequence files
# Author: Ryan D. Najac
# Date of modification: 2024-05-03

function multifastqc {
  local input_dir="$1"
  local output_dir="$2"

  mkdir -pv "$output_dir" || return 1
  fastqc -o "$output_dir" --noextract --memory 1024 -t "$(nproc)" "$input_dir"/*
  zip -r "${output_dir}/fastqc_html.zip" "$output_dir"/*.html
}

readonly PATH_TO_FILES="/home/ubuntu/rnaseq/"
readonly PROJECT_FOLDERS=("agx51" "ra")

for folder in "${PROJECT_FOLDERS[@]}"; do
  fastqc_output="${PATH_TO_FILES}/${folder}/fastqc"
  multifastqc "${PATH_TO_FILES}/${folder}/fastq" "$fastqc_output"
  echo "Finished processing $folder"
done

#======================= CHECKSUM FUNCTIONS =============================

# Generate checksums for all files in the current directory
generate_checksums() {
  local extension=${1:-"fastq.gz"}
  local checksum_file="md5sums.txt"
  echo "Generating checksums for *.$extension files in the current directory..."
  for file in *.$extension; do
    echo "$(md5sum "$file" | cut -d ' ' -f 1)\t$file" >> $checksum_file
  done
  echo "Checksums stored in $checksum_file"
}

# Verify checksums for all files in the current directory
verify_checksums() {
  local md5sums_txt="${1:-md5sums.txt}"
  echo "Verifying checksums..."
  while read -r expected_checksum file; do
    local status="[ OK ]"
    [[ "$(md5sum "$file" | cut -d ' ' -f 1)" != "$expected_checksum" ]] && status="[FAIL]"
    echo "${status} ${file}"
  done < $md5sums_txt
}
