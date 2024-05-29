#!/bin/bash
#
# USE WITH CAUTION
# Function to sync the folder containing this script with a remote server, ignoring hidden files and folders.

rsync -avz --progress --exclude='.*' "$(dirname "$(realpath "$0")")/" "aws:~"
