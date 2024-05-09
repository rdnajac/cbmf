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

update_ssh_config "$1"

