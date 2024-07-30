#!/bin/bash
#
## Keep this folder in sync with a remote machine
set -euo pipefail
set -x

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# This is set in ~/.ssh/config (`Host my-ec2`)
REMOTE_ALIAS="my-ec2"

echo $REMOTE_ALIAS

if [[ $(uname) == "Darwin" ]]; then
	# assume were on the host if we're on a mac
	# XXX: this is not a good assumption
	rsync -avz --delete "$SCRIPTS_DIR"/ "${REMOTE_ALIAS}:~/scripts/"

	# DO NOT delete files on the remote that are not on the local
	# (i.e. only copy files from the remote that are not on the local)
	# rsync -avz "${REMOTE_ALIAS}:~/out/" "$LOCAL_DATA_IN"/
else
	# change behavior if this script is run on the remote
	
	# add the scripts directory to the PATH if it's not already there
	if (grep -q "export PATH=\$PATH:~/scripts" ~/.bashrc); then
		echo "export PATH=\$PATH:~/scripts" >> ~/.bashrc
	fi

fi

