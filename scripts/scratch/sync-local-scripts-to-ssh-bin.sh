#!/bin/bash
#
# Use this file to sync local scripts to the remote server's bin directory
# This script will copy all files in the local bin directory to the remote server's bin directory

# use my-ec2, which is defined in ~/.ssh/config
REMOTE_SERVER=my-ec2
DEST_DIR=/usr/local/bin
DEST_PATH=${REMOTE_SERVER}:${DEST_DIR}

THIS_DIR=$(dirname $0)
# Copy all .sh files in the local bin directory to the remote server's bin directory
for file in $THIS_DIR/*.sh; do
    echo "Copying $file to $DEST_PATH"
    # Use rsync instead of scp for better performance and to preserve permissions
    rsync -e "ssh -o StrictHostKeyChecking=no" --rsync-path="sudo rsync" "$file" "$DEST_PATH"
done
