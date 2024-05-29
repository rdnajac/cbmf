#!/bin/bash
#
# USE WITH CAUTION
# Function to sync the folder containing this script with a remote server.

rsync -avz --progress "$(dirname "$(realpath "$0")")" "aws:~"
