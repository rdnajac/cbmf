
# bowtie2 -p "$threads" --mm -x "$ref" -1 "$R1" -2 "$R2" -S # samtools view -buS "$sam_file" | samtools sort -@ "$threads" -o "$sorted_bam"

# ``` sh
# bowtie2 -p 8 --mm -x index/lambda_virus \
#         -1 reads/reads_1.fq -2 reads/reads_2.fq \
#         -p 8 | samtools view -bS - | samtools sort -o sorted.bam
# ```


# This makefile should basically do all that
# we expect 1 of two suffices: `R{1..2}_001.fastq.gz` please use brace exapnsion 
# R1FASTQ := $(wildcard *_R1_001.fastq.gz)
# R2FASTQ := $(wildcard *_R2_001.fastq.gz)
THREADS := $(shell nproc)
REF ?= $(HUMANREF) 
# REF ?= $(MOUSEREF) 
SAMPLE_DIR ?= /path/to/directory
SAMPLES := $(basename $(notdir $(wildcard $(SAMPLE_DIR)/*_R1_001.fastq.gz)))

BOWTIE2 := bowtie2 -p $(THREADS) -mm -x $(REF) -1 $< -2 $(word 2, $^) -S $@
#process_samples "$1" "$2"

# given a sample name, do the bowtie2

Default target
all: $(SAMPLES:%=$(SAMPLE_DIR)/%.sorted.bam)

# Rule to perform alignment and generate SAM files
$(SAMPLE_DIR)/%.sam: $(SAMPLE_DIR)/%_R1_001.fastq.gz $(SAMPLE_DIR)/%_R2_001.fastq.gz
@echo "Aligning $*"
bowtie2 -p $(THREADS) --mm -x $(REF) -1 $< -2 $(word 2, $^) -S $@

Rule to perform alignment and generate sorted BAM files
$(SAMPLE_DIR)/%.sorted.bam: $(SAMPLE_DIR)/%.sam
@echo "Converting SAM to sorted BAM for $*"
samtools view -buS $
