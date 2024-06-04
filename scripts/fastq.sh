#!/bin/bash

# This script handles moving the fastq data around
#
# Azenta provides checksum files for each fastq. check them after download
for f in *.fastq.gz; do
  md5sum -c "${f}.md5"
done

# again but one line
for f in *.fastq.gz; do md5sum -c "${f}.md5"; done

# again, bit with a condition that if its all good we self destrict and delete all the md5 files
for f in *.fastq.gz; do
  md5sum -c "${f}.md5" && rm "${f}.md5" || echo "Checksum failed for ${f}"
done

# one liner
for f in *.fastq.gz; do md5sum -c "${f}.md5" && rm -f "${f}.md5" || echo "Checksum failed for ${f}"; done

# create a single checksum file for all the fastq files
md5sum ./*.fastq.gz > checksums.md5

# print it
cat checksums.md5
