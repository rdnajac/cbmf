#!/bin/bash
# bash functions

sync_fastq() {
    # Sync all fastq.gz files in the current folder to the S3 URI provided as the first argument
    if [[ -z "$1" ]]; then
        echo "Usage: sync_fastq <s3-uri>"
        return 1
    fi

    aws s3 sync . "$1" --exclude "*" --include "*.fastq.gz" --exclude "*/"
}
