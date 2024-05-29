#!/bin/bash

# Variables - replace these with your actual file paths or URLs
REFERENCE_GENOME_URL="ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REFERENCE_GENOME_FA="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE_URL="ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"
GTF_FILE="Homo_sapiens.GRCh38.84.gtf"
SNP_FILE_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz"
SNP_FILE="snp144Common.txt"

# Download and prepare the reference genome sequence
wget $REFERENCE_GENOME_URL
gzip -d ${REFERENCE_GENOME_FA}.gz
mv $REFERENCE_GENOME_FA genome.fa

# Download and prepare the GTF file
wget $GTF_FILE_URL
gzip -d ${GTF_FILE}.gz
mv $GTF_FILE genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon

# Download and prepare the SNP file
wget $SNP_FILE_URL
gzip -d ${SNP_FILE}.gz
awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' $SNP_FILE > ${SNP_FILE}.ensembl
hisat2_extract_snps_haplotypes_UCSC.py genome.fa ${SNP_FILE}.ensembl genome

# Build HFM index
hisat2-build -p 16 genome.fa genome

# Build HGFM index with SNPs
hisat2-build -p 16 --snp genome.snp --haplotype genome.haplotype genome.fa genome_snp

# Build HGFM index with transcripts
hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa genome_tran

# Build HGFM index with SNPs and transcripts
hisat2-build -p 16 --snp genome.snp --haplotype genome.haplotype --exon genome.exon --ss genome.ss genome.fa genome_snp_tran
