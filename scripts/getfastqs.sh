#!/bin/bash
BUCKET='lab-aaf-ngs-data-archive'
folder1="ChIPseq/20230216_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_CB/00_fastq/"
folder2="ChIPseq/20230308_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_rep2redo_CB/00_fastq/"
folder3="ChIPseq/20230308_ChIP-seq_HEL_STAT5_ruxolitinib_SGC-CBP30_CB/00_fastq/"

folders=("$folder1" "$folder2" "$folder3")

# # Retrieve a list of all objects in Glacier or Glacier Deep Archive storage class within the specified path
# FOLDER3OBJECTS=$(aws s3api list-objects-v2 --bucket "$BUCKET" --prefix "$folder3" --query "Contents[?StorageClass=='GLACIER' || StorageClass=='DEEP_ARCHIVE'].Key" --output text)
# restore the objects
# for OBJECT in $FOLDER3OBJECTS; do
# aws s3api head-object --bucket "$BUCKET" --key "$OBJECT" --query "Restore" --output text
# done

# files saved in ./run1 ./run2 and ./run3
