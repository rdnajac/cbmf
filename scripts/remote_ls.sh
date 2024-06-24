#!/bin/bash

# Global variable for remote host alias
remoteHostAlias="my-ec2"
scriptName=$(basename "$0")

# Function to run the script remotely and broadcast the output to all terminals
remote_execute()
{
	local output

	# Copy the script to the remote server
	scp "$scriptName" "$remoteHostAlias:~/$scriptName"

	# Run the script on the remote server and capture the output
	output=$(ssh "$remoteHostAlias" "~/$scriptName")

	# Broadcast the output to all terminals on the remote server
	ssh -t -t "$remoteHostAlias" "sudo wall <<< '$output'"
}

# Usage example
remote_execute
