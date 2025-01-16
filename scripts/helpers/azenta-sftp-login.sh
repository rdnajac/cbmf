#!/bin/bash
## Download data from Azenta via SFTP

# replace with your username
USER="hm2882_cumc_columbia"
# omit the `sftp://` prefix
HOST="sftp.genewiz.com"
PORT="22"

sftp -P "$PORT" "$USER@$HOST"
# sftp -P 22 hm2882_cumc_columbia@sftp.genewiz.com

# enter `yes` to accept the key
# enter the password
# files will be downloaded to the current directory from where the script is run
# use `lcd` to change the local directory (download location)
# run `get -fpR *` to download all files and directories
# this might take a while
