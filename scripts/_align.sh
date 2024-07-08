#!/bin/bash
#
## Wrapper script for paired-end read alignment using HISAT2 or Bowtie2.

# #######################################
# Runs either bowtie2 or hisat2 on a pair of fastq files
# Usage:
#  $ _align.sh <aligner> <r1_fastq> <r2_fastq> <reference>
# Globals:
#   None
# Arguments:
#  $1: /path/to/aligner (hisat2 or bowtie2)
#  $2: read 1 fastq file
#  $3: read 2 fastq file
#  $4: reference genome
#
# Outputs:
#   STDOUT: sam file
#   STDERR: log file
#
# Notes:
#   Callee must handle capturing output!
#   Uses -p to specify the number of threads to use.
#   Uses --mm use memory-mapped I/O for performance.
#######################################
main() {
	_align "$@"
	"$(exit $?)" || echo "Alignment failed"
}

main "$@"
