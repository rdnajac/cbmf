#!/bin/bash
#
# Copies a shell script to my-ec2:~/ and executes it remotely.

# Usage: ./remote_exec.sh <script.sh>
# Example: ./remote_exec.sh my_script.sh

# Check if script is provided
# If not, print usage and exit
if [ $# -eq 0 ]; then
	echo "Usage: $0 <script.sh>"
	exit 1
fi

# Check if script exists
# If not, print error and exit
if [ ! -f "$1" ]; then
	echo "Error: $1 does not exist"
	exit 1
fi

# Copy script to my-ec2:~/
scp "$1" my-ec2:~/

# Execute script remotely
ssh my-ec2 "bash ~/$1"

# Done
exit 0
