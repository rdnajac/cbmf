#!/bin/bash
## Given a directory, consolidate all the md5 files into a single file
set -euo pipefail
set -x

# Check if the user has provided a directory
if [ $# -ne 1 ]; then
	echo "Usage: $0 <directory>"
	exit 1
fi

cd "$1" || {
	echo "Directory does not exist"
	exit 2
}

# if there is no file called "md5sums.txt", create one
if [ ! -f md5sums.txt ]; then
	cat ./*.md5 > md5sums.txt
fi
# check the md5sum of the consolidated file
md5sum -c md5sums.txt

# # ask user if they want to delete the individual md5 files
# read -r -p "Do you want to delete the individual md5 files? [y/N] " response

# if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
# 	rm -v ./*.md5
# fi
