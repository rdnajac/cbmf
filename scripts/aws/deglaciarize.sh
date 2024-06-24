#!/bin/bash

# This file is used to de-glaciarize files in a specific S3 bucket.

# define variables 
BUCKET="lab-aaf-ngs-data-archive"

FOLDER="RNAseq/20220807_DNMT3A_RHOA_PreLymphoma_APL/"
FILES="nohup.out runSTAR.sh runSTAR_v2.sh runfeatureCounts.sh runfeatureCounts2.sh samples.txt"
KEYS=()
for FILE in $FILES; do
  KEYS+=("$FOLDER$FILE")
done

# define the length of time to de-glaciarize the files for
DAYS=7

restore() {
  aws s3api restore-object --bucket $BUCKET --key $1 --restore-request '{"Days":'$DAYS',"GlacierJobParameters":{"Tier":"Standard"}}' | tee -a deglaciarize.log
}

check() {
  aws s3api head-object --bucket $BUCKET --key $1 | grep Restore | awk -v key="$1" '{print key, $0}' | tee -a deglaciarize.log
}

for KEY in "${KEYS[@]}"; do
  # restore $KEY
  check $KEY
done 

# one-liner to aws sync that to a local folder:
# aws s3 sync --force-glacier-transfer s3://lab-aaf-ngs-data-archive/ /path/to/local/folder/
