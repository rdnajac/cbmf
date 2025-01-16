#!/bin/bash
# syncs all fastq.gz files in current folder to s3URI $1
# aws s3 sync /path/to/source s3://your-bucket-name/your-prefix/ --exclude "*" --include "*.fastq.gz" --exclude "*/"

aws s3 sync . "$1" --exclude "*" --include "*.fastq.gz" --exclude "*/"
