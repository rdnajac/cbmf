#!/bin/bash
#
# RNAseq analysis pipeline

FILE_PATH=$(realpath "$0")
FILE_NAME=$(basename "$0")

REMOTE_HOST=my-ec2
EXECUTED=false

remote_exec()
{
	scp "$FILE_PATH" "$REMOTE_HOST:~/" && ssh "$REMOTE_HOST" "bash ~/${FILE_NAME}"
}

# THE COMMANDS BELOW WILL BE EXECUTED ON THE REMOTE HOST
cd /FASTQ || exit
# aws s3 sync s3://lab-aaf-ngs-data-archive/RNAseq/20231127_Tet2Rhoa-sgID2-3-cloneB_APL/ . &
# aws s3 sync s3://lab-aaf-ngs-data-archive/RNAseq/20240409_Tet2Rhoa-S1P1_RA/ . &
# aws s3 sync s3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2Rhoa-AGX51_APL/ . &
touch hi

# Call remote_exec only once
if [ "$EXECUTED" = false ]; then
	remote_exec
	EXECUTED=true
fi
