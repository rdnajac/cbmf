#!/bin/bash

## NAME:
##   cbmf - Command-line Bioinformatics Multi-Function tool
##
## SYNOPSIS:
##   cbmf <command> [options] [arguments]
##
## DESCRIPTION:
##   This script provides various bioinformatics functions including alignment,
##   quality control, genome initialization, and file format conversion.
##
## LAST REVISION: 2024-08-05

set -o errexit
set -o nounset

# Global variables
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
VALID_SPECIES=("human" "mouse")
VALID_ALIGNERS=("bwa" "hisat2" "bowtie2" "star" "subread-align" "subjunct")
GENOMES_DIR="$HOME/genomes"
GENOMES_MIRROR="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"
WGET_OPTS="--no-verbose --no-parent --no-directories"
MOUSE_GCA="000/001/635/GCA_000001635.9_GRCm39"
HUMAN_GCA="000/001/405/GCA_000001405.15_GRCh38"
MOUSE_PREFIX="GCA_000001635.9_GRCm39_full_analysis_set"
HUMAN_PREFIX="GCA_000001405.15_GRCh38_full_analysis_set"

# Helper functions
pr_info() {
    printf "\e[1;34m[INFO] %s\e[0m\n" "$1"
}

pr_error() {
    printf "\e[1;31m[ERROR] %s\e[0m\n" "$1" >&2
}

bail() {
    pr_error "$1"
    exit 1
}

validate() {
    local value="$1"
    local -n valid_options="$2"
    value=$(tr '[:upper:]' '[:lower:]' <<< "$value")
    for option in "${valid_options[@]}"; do
        if [ "$value" = "$option" ]; then
            return 0
        fi
    done
    return 1
}

# Alignment functions
align_reads() {
    local aligner="$1"
    local r1_fastq="$2"
    local r2_fastq="$3"
    local reference="$4"

    case "$aligner" in
        hisat2)
            hisat2 -p "$(nproc)" --mm -1 "$r1_fastq" -2 "$r2_fastq" -x "$reference" --dta
            ;;
        bowtie2)
            bowtie2 -p "$(nproc)" --mm -1 "$r1_fastq" -2 "$r2_fastq" -x "$reference"
            ;;
        *)
            bail "Invalid aligner: $aligner"
            ;;
    esac
}

align_to_bam() {
    local aligner="$1"
    local r1_fastq="$2"
    local r2_fastq="$3"
    local reference="$4"
    local bam_file="${r1_fastq%_R1*}_sorted_markdup.bam"

    align_reads "$aligner" "$r1_fastq" "$r2_fastq" "$reference" \
        | samtools sort -n -@ "$(nproc)" - \
        | samtools fixmate -@ "$(nproc)" -m - - \
        | samtools sort -@ "$(nproc)" - \
        | samtools markdup -@ "$(nproc)" - "$bam_file"
}

# Genome initialization functions
download_index() {
    local species="$1"
    local prefix="$2"
    local gca="$3"
    local index="$4"

    local file="${prefix}.fna.${index}_index.tar.gz"
    wget $WGET_OPTS -O "$file" "$GENOMES_MIRROR/$gca/seqs_for_alignment_pipelines.ucsc_ids/$file"
    tar -xzvf "$file"
    rm "$file"
    pr_info "${index} index for $species downloaded and extracted."
}

download_fasta() {
    local species="$1"
    local prefix="$2"
    local gca="$3"

    wget $WGET_OPTS "$GENOMES_MIRROR/$gca/seqs_for_alignment_pipelines.ucsc_ids/${prefix}.fna.fai"
    wget $WGET_OPTS -O "${prefix}.fna.gz" "$GENOMES_MIRROR/$gca/seqs_for_alignment_pipelines.ucsc_ids/${prefix}.fna.gz"
    gunzip "${prefix}.fna.gz"
    pr_info "FASTA file for $species downloaded and extracted."
}

download_refseq() {
    local species="$1"
    local prefix="$2"
    local gca="$3"
    local format="$4"

    wget $WGET_OPTS -O "${prefix}.refseq_annotation.${format}.gz" "$GENOMES_MIRROR/$gca/seqs_for_alignment_pipelines.ucsc_ids/${prefix}.refseq_annotation.${format}.gz"
    gunzip "${prefix}.refseq_annotation.${format}.gz"
    pr_info "RefSeq annotation (${format}) for $species downloaded and extracted."
}

init_genome() {
    local species="$1"
    shift
    local prefix
    local gca

    case "$species" in
        mouse)
            prefix="$MOUSE_PREFIX"
            gca="$MOUSE_GCA"
            ;;
        human)
            prefix="$HUMAN_PREFIX"
            gca="$HUMAN_GCA"
            ;;
        *)
            bail "Invalid species. Use 'mouse' or 'human'."
            ;;
    esac

    mkdir -p "$GENOMES_DIR/$species"
    cd "$GENOMES_DIR/$species" || bail "Failed to change directory"

    if [ $# -eq 0 ]; then
        download_fasta "$species" "$prefix" "$gca"
        download_index "$species" "$prefix" "$gca" "bowtie"
        download_index "$species" "$prefix" "$gca" "bwa"
        download_index "$species" "$prefix" "$gca" "hisat2"
        download_refseq "$species" "$prefix" "$gca" "gff"
        download_refseq "$species" "$prefix" "$gca" "gtf"
    else
        for component in "$@"; do
            case "$component" in
                fasta)
                    download_fasta "$species" "$prefix" "$gca"
                    ;;
                bowtie|bwa|hisat2)
                    download_index "$species" "$prefix" "$gca" "$component"
                    ;;
                gff|gtf)
                    download_refseq "$species" "$prefix" "$gca" "$component"
                    ;;
                *)
                    pr_error "Invalid component: $component. Skipping."
                    ;;
            esac
        done
    fi

    pr_info "Requested $species genome components downloaded and extracted successfully."
}

# File conversion functions
bam_to_cram() {
    local bam_file="$1"
    local id="$(basename "$bam_file" .bam)"
    local cram_file="${id}.cram"
    local reference="/home/ubuntu/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set.fna.bgz"

    samtools view -@"$(nproc)" --cram -T "$reference" "$bam_file" > "$cram_file"
    pr_info "Converted $bam_file to $cram_file"
}

cram_to_bam() {
    local cram_file="$1"
    local id="$(basename "$cram_file" .cram)"
    local bam_file="${id}.bam"

    if [ -f "$bam_file" ]; then
        mv "$bam_file" "${bam_file}.bak"
    fi

    samtools view -@"$(nproc)" --bam "$cram_file" > "$bam_file"
    pr_info "Converted $cram_file to $bam_file"
}

# Quality control function
run_fastqc() {
    local out_dir="${1:-fastqc}"
    local memsize=4096

    mkdir -vp "$out_dir"
    fastqc -t "$(nproc)" --memory "$memsize" --noextract -o "$out_dir" ./*
    pr_info "FastQC analysis completed. Results in $out_dir"
}

# Main command functions
cmd_align() {
    local species=""
    local aligner=""
    local input_dir=""
    local output_dir=""

    while getopts ":s:a:i:o:h" opt; do
        case $opt in
            s) species="$OPTARG" ;;
            a) aligner="$OPTARG" ;;
            i) input_dir="$OPTARG" ;;
            o) output_dir="$OPTARG" ;;
            h) cmd_align_help; exit 0 ;;
            \?) bail "Invalid option: -$OPTARG" ;;
        esac
    done

    [ -z "$species" ] && bail "Species (-s) is required for alignment"
    [ -z "$aligner" ] && bail "Aligner (-a) is required for alignment"
    [ -z "$input_dir" ] && bail "Input directory (-i) is required"
    [ -z "$output_dir" ] && bail "Output directory (-o) is required"

    validate "$species" VALID_SPECIES || bail "Invalid species: $species"
    validate "$aligner" VALID_ALIGNERS || bail "Invalid aligner: $aligner"

    local reference="$GENOMES_DIR/${species}/${species^}_genome"
    mkdir -p "$output_dir"

    for r1_fastq in "$input_dir"/*_R1.fastq.gz; do
        r2_fastq="${r1_fastq/_R1/_R2}"
        [ -f "$r2_fastq" ] || bail "Paired read file not found: $r2_fastq"
        align_to_bam "$aligner" "$r1_fastq" "$r2_fastq" "$reference"
    done

    pr_info "Alignment completed. BAM files are in $output_dir"
}

cmd_init() {
    [ $# -lt 1 ] && bail "Usage: cbmf init {mouse|human} [components...]"
    init_genome "$@"
}

cmd_qc() {
    local input_dir="${1:-.}"
    run_fastqc "$input_dir"
}

cmd_convert() {
    [ $# -ne 1 ] && bail "Usage: cbmf convert <file.bam|file.cram>"
    local file="$1"
    case "$file" in
        *.bam) bam_to_cram "$file" ;;
        *.cram) cram_to_bam "$file" ;;
        *) bail "Unknown file type: $file" ;;
    esac
}

cmd_help() {
    cat << EOF
Usage: cbmf <command> [options] [arguments]

Commands:
  align         Perform sequence alignment
  init          Initialize genome files
  qc            Perform quality control
  convert       Convert between BAM and CRAM formats
  help          Show this help message

Run 'cbmf help <command>' for more information on a specific command.
EOF
}

cmd_align_help() {
    cat << EOF
Usage: cbmf align [options]

Perform sequence alignment.

Options:
  -s SPECIES    Specify species (human or mouse)
  -a ALIGNER    Specify aligner (bwa, hisat2, bowtie2, star, subread-align, subjunct)
  -i INPUT      Input directory containing FASTQ files
  -o OUTPUT     Output directory for alignment results
  -h            Show this help message
EOF
}

# Main execution
main() {
    [ $# -eq 0 ] && { cmd_help; exit 1; }

    local cmd="$1"
    shift

    case "$cmd" in
        align) cmd_align "$@" ;;
        init) cmd_init "$@" ;;
        qc) cmd_qc "$@" ;;
        convert) cmd_convert "$@" ;;
        help)
            if [ $# -eq 0 ]; then
                cmd_help
            else
                case "$1" in
                    align) cmd_align_help ;;
                    init|qc|convert) pr_info "$1 help not yet implemented" ;;
                    *) bail "Unknown command: $1" ;;
                esac
            fi
            ;;
        *) bail "Unknown command: $cmd" ;;
    esac
}

main "$@"
