#!/bin/bash

update_ssh_config() {
    new_hostname="$1"
    config_file="$HOME/.ssh/config"

    if [ -z "$new_hostname" ]; then
        echo "Usage: update_ssh_config new_hostname"
        return 1
    fi

    if [ ! -f "$config_file" ]; then
        echo "SSH config file does not exist: $config_file"
        return 1
    fi

    # Modify the sed command to work on macOS by adding an empty string after -i
    sed -i '' "\$s/^  HostName .*/  HostName $new_hostname/" "$config_file"

    echo "Hostname in $config_file has been updated to $new_hostname."
}

# Assuming the instance ID is stored in the MY_AWS_ID environment variable...
update_ssh_config "$(aws ec4 describe-instances --instance-ids "${MY_AWS_ID}" --query 'Reservations[0].Instances[0].PublicDnsName' --output text)"

# TODO: check if the instance is running, start it if it is not

# if its nor already running:
# update_ssh_config "$(aws ec2 start-instances --instance-ids "${MY_AWS_ID}" --query 'StartingInstances[0].InstanceId' --output text)"
