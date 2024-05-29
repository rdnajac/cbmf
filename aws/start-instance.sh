#!/bin/bash

# Start the instance and capture the instance ID
#instance_id=$(aws ec2 start-instances --instance-ids <your_instance_id> --query 'StartingInstances[0].InstanceId' --output text)
instance_id="${MY_AWS_ID}"

# Wait for the instance to be in the running state
aws ec2 wait instance-running --instance-ids "$instance_id"

# Retrieve the public DNS of the running instance
public_dns=$(aws ec2 describe-instances --instance-ids "$instance_id" --query 'Reservations[0].Instances[0].PublicDnsName' --output text)

echo "Public DNS: $public_dns"
# if update-aws-ssh.sh exists and is executable, run it
[[ -x update-aws-ssh.sh ]] && ./update-aws-ssh.sh "$public_dns" && echo "SSH updated"
