#!/bin/bash
## Upload data to AWS S3 bucket
## This script assumes you have the AWS CLI installed and configured

# Set the bucket name
PROJECT="YYYYMMDD-PROJECT-UNI"
S3_BUCKET="s3://PATH/TO/BUCKET/$PROJECT"

# Azenta project layout (after sftp download)
# .
# ├── 30-##########
# │   ├── 00_fastq/
# │   ├── *.fastq.gz
# │   └── *.fastq.gz.md5
# └── Azenta_30-##########_Data_Report.html

# Change as needed
AZENTAID="30-##########"
FASTQ_DIR="00_fastq"


