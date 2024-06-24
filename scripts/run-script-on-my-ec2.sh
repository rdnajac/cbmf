#!/bin/bash
#
# Copies a shell script to my-ec2:~/ and executes it remotely.

# Given a file...
if [ "$1" = "" ]; then
	echo "Usage: $0 <script>"
	exit 1
fi

# execute it over ssh on my-ec2
ssh my-ec2 "bash -s" <"$1" &
# do I need -t? no. -t -t? I don't think so. -T?
#
