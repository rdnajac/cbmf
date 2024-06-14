#!/bin/bash
#
# This script is a pipeline to align paired-end reads to a reference genome

MY_DIR=$(dirname "${0}")

# Source the configuration file
source "${MY_DIR}/utils.sh"

# Functions ______________________________________________________________

fastq2bam() {
  local id="$1"
  local suffix="$2"  # including .fastq.gz
  local r1="${id}_R1${suffix}"
  local r2="${id}_R2${suffix}"
  local log_file="${id}_hisat2.log"
  local bam_file="${id}_sorted_markdup.bam"

  # Optionally, replace the first line of the pipeline
  # bowtie2 -mm --threads "$THREADS" -x "$REF" -1 "$r1" -2 "$r2" 2>> "${id}.log" \
  
  hisat2 --mm --threads "$THREADS" -x "$HISAT_REF_INDEX" -1 "$r1" -2 "$r2" --met-stderr 2> >(tee "$log_file" >&2) \
  | samtools sort    -@ "$THREADS" -n -   \
  | samtools fixmate -@ "$THREADS" -m - - \
  | samtools sort    -@ "$THREADS"    -   \
  | samtools markdup -@ "$THREADS"    - "$bam_file"


  # write protect all *.fastq.gz and *.bam files
  chmod a-w ./*.bam
}

MERGED_GTF="stringtie_merged.gtf"


# StringTie wrapper
# @param $1 bam file
# @param $2 reference (.gtf or .gff)
run_stringtie() {
  local sample="$1"
  local output="${sample%.bam}.gtf"
  local reference="$2"

  $STRINGTIE -p "$THREADS" -G "$reference" -o "$output" "$sample"
}

# Given Run strintie on every bamfile in this firectory
do_run_stringtie() {
  for f in ./*.bam; do
    run_stringtie "$f" "$REFSEQ_ANNOTATION_GTF" &
  done
  wait
  echo "All samples processed with StringTie."
}

# Function to merge GTF files
# operates in the current directory
merge_gtfs() {
    stringtie --merge -p "$THREADS" -G "$REFSEQ_ANNOTATION_GTF" -o "$MERGED_GTF" ./*.gtf 
    echo "GTF files merged into $MERGED_GTF."
}

# Function to process samples with merged GTF
process_with_merged_gtf() {
    for f in ./*.bam; do
      local id="${f%.bam}"
      local folder="ballgown_input/$id"
      mkdir -vp "$folder"
      stringtie -p "$THREADS" -G "$MERGED_GTF" -e -b "$folder" -o "${id}_stringtie.gtf" "$f" &
    done
    wait
    echo "All samples processed with StringTie using merged GTF."
}

# assume were in the folder with the bam files ready for stringtie processsing
do_run_stringtie
merge_gtfs
process_with_merged_gtf
echo "All done."
exit 0


# Main logic
# run_stringtie
# merge_gtfs
# process_with_merged_gtf

# StringTie v2.2.3 usage:
# stringtie <in.bam ..> [-G <guide_gff>] [-l <prefix>] [-o <out.gtf>] [-p <cpus>]
#  [-v] [-a <min_anchor_len>] [-m <min_len>] [-j <min_anchor_cov>] [-f <min_iso>]
#  [-c <min_bundle_cov>] [-g <bdist>] [-u] [-L] [-e] [--viral] [-E <err_margin>]
#  [--ptf <f_tab>] [-x <seqid,..>] [-A <gene_abund.out>] [-h] {-B|-b <dir_path>}
#  [--mix] [--conservative] [--rf] [--fr]
# Assemble RNA-Seq alignments into potential transcripts.
# Options:
#  --version : print just the version at stdout and exit
#  --conservative : conservative transcript assembly, same as -t -c 1.5 -f 0.05
#  --mix : both short and long read data alignments are provided
#         (long read alignments must be the 2nd BAM/CRAM input file)
#  --rf : assume stranded library fr-firststrand
#  --fr : assume stranded library fr-secondstrand
#  -G reference annotation to use for guiding the assembly process (GTF/GFF)
#  --ptf : load point-features from a given 4 column feature file <f_tab>
#  -o output path/file name for the assembled transcripts GTF (default: stdout)
#  -l name prefix for output transcripts (default: STRG)
#  -f minimum isoform fraction (default: 0.01)
#  -L long reads processing; also enforces -s 1.5 -g 0 (default:false)
#  -R if long reads are provided, just clean and collapse the reads but
#     do not assemble
#  -m minimum assembled transcript length (default: 200)
#  -a minimum anchor length for junctions (default: 10)
#  -j minimum junction coverage (default: 1)
#  -t disable trimming of predicted transcripts based on coverage
#     (default: coverage trimming is enabled)
#  -c minimum reads per bp coverage to consider for multi-exon transcript
#     (default: 1)
#  -s minimum reads per bp coverage to consider for single-exon transcript
#     (default: 4.75)
#  -v verbose (log bundle processing details)
#  -g maximum gap allowed between read mappings (default: 50)
#  -M fraction of bundle allowed to be covered by multi-hit reads (default:1)
#  -p number of threads (CPUs) to use (default: 1)
#  -A gene abundance estimation output file
#  -E define window around possibly erroneous splice sites from long reads to
#     look out for correct splice sites (default: 25)
#  -B enable output of Ballgown table files which will be created in the
#     same directory as the output GTF (requires -G, -o recommended)
#  -b enable output of Ballgown table files but these files will be
#     created under the directory path given as <dir_path>
#  -e only estimate the abundance of given reference transcripts (requires -G)
#  --viral : only relevant for long reads from viral data where splice sites
#     do not follow consensus (default:false)
#  -x do not assemble any transcripts on the given reference sequence(s)
#  -u no multi-mapping correction (default: correction enabled)
#  -h print this usage message and exit
#  --ref/--cram-ref reference genome FASTA file for CRAM input

# Transcript merge usage mode:
#   stringtie --merge [Options] { gtf_list | strg1.gtf ...}
# With this option StringTie will assemble transcripts from multiple
# input files generating a unified non-redundant set of isoforms. In this mode
# the following options are available:
#   -G <guide_gff>   reference annotation to include in the merging (GTF/GFF3)
#   -o <out_gtf>     output file name for the merged transcripts GTF
#                     (default: stdout)
#   -m <min_len>     minimum input transcript length to include in the merge
#                     (default: 50)
#   -c <min_cov>     minimum input transcript coverage to include in the merge
#                     (default: 0)
#   -F <min_fpkm>    minimum input transcript FPKM to include in the merge
#                     (default: 1.0)
#   -T <min_tpm>     minimum input transcript TPM to include in the merge
#                     (default: 1.0)
#   -f <min_iso>     minimum isoform fraction (default: 0.01)
#   -g <gap_len>     gap between transcripts to merge together (default: 250)
#   -i               keep merged transcripts with retained introns; by default
#                    these are not kept unless there is strong evidence for them
#   -l <label>       name prefix for output transcripts (default: MSTRG)
# ubuntu@ip-172-31-32-180:~/apl$

