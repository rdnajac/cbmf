#!/bin/bash
## Download data from Azenta via SFTP

# Set the variables
USER="hm2882_cumc_columbia"
# replace with your username
HOST="sftp.genewiz.com"
# omit the `sftp://` prefix
PORT="22"

sftp -P "$PORT" "$USER@$HOST"
# enter `yes` to accept the key
# enter the password
# files will be downloaded to the current directory from where the script is run
# use `lcd` to change the local directory (download location)
# run `get -fpR *` to download all files and directories
# this might take a while depending on the size of the data
