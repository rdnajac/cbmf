#!/bin/bash

# This script checks the environment variable MY_EC2_ID for an AWS instance ID.
# It starts the instance if it is not already running and updates the SSH
# configuration file with the instance's public DNS name.

# Exit if MY_EC2_ID is not set
[[ -z "${MY_EC2_ID}" ]] && echo "MY_EC2_ID is not set." >&2 && exit 1

# Check if the instance is running, start it if it is not
if [ "$(aws ec2 describe-instances --instance-ids "${MY_EC2_ID}" --query 'Reservations[0].Instances[0].State.Name' --output text)" != "running" ]; then
  echo "Starting the instance..."
  aws ec2 start-instances --instance-ids "${MY_EC2_ID}"
fi

# Get the public DNS name of the instance
new_hostname="$(aws ec2 describe-instances --instance-ids "${MY_EC2_ID}" --query 'Reservations[0].Instances[0].PublicDnsName' --output text)"

# Define the SSH config file path
config_file="$HOME/.ssh/config"

# Check if the SSH config file exists
if [ ! -f "$config_file" ]; then
  echo "SSH config file does not exist: $config_file" >&2
  exit 2
fi

# Update the SSH config file with the new hostname
sed -i '' "\$s/^  HostName .*/  HostName $new_hostname/" "$config_file"

echo "Hostname in $config_file has been updated to $new_hostname."

tail -n 4 $config_file

