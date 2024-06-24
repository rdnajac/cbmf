#!/bin/bash
#
# Assemble transcripts with stringtie
# https://github.com/gpertea/stringtie
#
# Input: folder with .bam files
#
# Workflow:
#     0a. Map the reads for each sample to the reference genome
#     0b. Sort and convert the SAM files to BAM
#     1. Assemble and quantify expressed genes and transcripts (StringTie)
#     2. Merge the GTF files from all samples (stringtie --merge)
#     3. Examine how the transcripts compare with the reference annotation (gffcompare)
#     4. Estimate transcript abundances and create table counts for Ballgown

set -euo pipefail
set -x

MY_DIR=$(dirname "${0}")
source "${MY_DIR}/utils.sh"
source "${MY_DIR}/stringtie.sh"



# Examine how the transcripts compare with the reference annotation (optional):$
# gffcompare –r chrX_data/genes/chrX.gtf –G –o merged stringtie_merged.gtf
# @brief Compare the transcripts with the reference annotation


do_run_stringtie "$MREF_GTF"
