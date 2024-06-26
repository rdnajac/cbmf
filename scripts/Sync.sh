#!/bin/bash
#
## These commands sync the scripts folder on my local machine
## to a folder on the remote machine (where that folder is in the PATH)
set -euo pipefail
set -x

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# remddote alias set in .ssh/config
REMOTE_ALIAS="my-ec2"

# define directories to keep in sync
LOCAL_DATA_DIR="${HOME}/cbmf/data"

# if ewre on darwin
if [[ $(uname) == "Darwin" ]]; then
	# assume were on the host
	# yes, delete files on the remote that are not on the local
	# (i.e. only copy files from the local that are not on the remote)
	rsync -avz --delete "$SCRIPTS_DIR"/ "${REMOTE_ALIAS}:~/scripts/"

	# DO NOT delete files on the remote that are not on the local
	# (i.e. only copy files from the remote that are not on the local)
	rsync -avz "${REMOTE_ALIAS}:~/out/" "$LOCAL_DATA_DIR"/
else
	if (grep -q "export PATH=\$PATH:~/scripts" ~/.bashrc); then
		echo "export PATH=\$PATH:~/scripts" >> ~/.bashrc
	fi

fi

# run a script over ssh
# runs the contents of a local file on a remote machine
ex-ec2() {
	ssh my-ec2 "bash -s" < "$1" &
}
