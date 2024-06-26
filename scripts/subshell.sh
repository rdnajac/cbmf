#!/bin/bash
#
## Wrapper to run commands in a subshell with notifications and logging

# function to send msg to sudo wall
broadcast()
{
	echo "$1" | sudo wall
}

# run commands in a subshell and redirect output from stderr to sudo wall
(
	# redirect stderr to stdout
	exec 2>&1

	# run the command
	"$@"

	# send notification
	broadcast "Command '$@' has finished"
) | sudo tee -a /var/log/subshell.log
