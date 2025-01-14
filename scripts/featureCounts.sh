#!/bin/bash
# JURKAT-CTRL-Clone1-LQ_sorted_markdup.bam
# JURKAT-CTRL-Clone2-LQ_sorted_markdup.bam
# JURKAT-CTRL-Clone3-LQ_sorted_markdup.bam
# JURKAT-PHF6-KO-Clone3-LQ_sorted_markdup.bam
# JURKAT-PHF6-KO-Clone7-LQ_sorted_markdup.bam
# JURKAT-PHF6-KO-Clone8-LQ_sorted_markdup.bam
# JURKAT-PHIP-KO-Clone10-MB_sorted_markdup.bam
# JURKAT-PHIP-KO-Clone22-LQ_sorted_markdup.bam
# JURKAT-PHIP-KO-Clone5-LQ_sorted_markdup.bam
# JURKAT-PHIP-KO-Clone9-MB_sorted_markdup.bam

ANNOTATION="/opt/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz"

FLAGS="-a ${ANNOTATION} -p --countReadPairs -T$(nproc)"

featureCounts "$FLAGS" -o total_counts.txt ./*.bam
featureCounts "$FLAGS" -o ctrl-phf6_counts.txt ./*CTRL*.bam ./*PHF6*.bam
featureCounts "$FLAGS" -o ctrl-phip_counts.txt ./*CTRL*.bam ./*PHIP*.bam
featureCounts "$FLAGS" -o phf6-phip_counts.txt ./*PHF6*.bam ./*PHIP*.bam
