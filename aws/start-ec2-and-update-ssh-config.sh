#!/bin/bash

# This script checks the environment variable MY_EC2_ID for an AWS instance ID.
# It starts the instance if it is not already running and updates the SSH
# configuration file with the instance's public DNS name.

# Exit if any command fails
set -euo pipefail

# Ensure MY_EC2_ID is set
: "${MY_EC2_ID:?MY_EC2_ID is not set.}"

# Get the instance state
instance_state=$(aws ec2 describe-instances --instance-ids "${MY_EC2_ID}" --query 'Reservations[0].Instances[0].State.Name' --output text)

# Start the instance if it is not running
if [[ "$instance_state" != "running" ]]; then
  echo "Starting the instance..."
  aws ec2 start-instances --instance-ids "${MY_EC2_ID}"
  echo "Waiting for the instance to start..."
  aws ec2 wait instance-running --instance-ids "${MY_EC2_ID}"
fi

# Get the public DNS name of the instance
new_hostname=$(aws ec2 describe-instances --instance-ids "${MY_EC2_ID}" --query 'Reservations[0].Instances[0].PublicDnsName' --output text)

# Define the SSH config file path
config_file="$HOME/.ssh/config"

# Ensure the SSH config file exists
if [[ ! -f "$config_file" ]]; then
  echo "SSH config file does not exist: $config_file" >&2
  exit 2
fi

# Update the SSH config file with the new hostname
sed -i '' "\$s/^  HostName .*/  HostName $new_hostname/" "$config_file"
#sed -i "\$s/^  HostName .*/  HostName $new_hostname/" "$config_file"

echo "Hostname in $config_file has been updated to $new_hostname."

# Display the last 4 lines of the SSH config file
tail -n 4 "$config_file"
