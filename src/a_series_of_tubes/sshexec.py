#!/usr/bin/env python3

import paramiko
import os
import logging


# Terminal colors
class bcolors:
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    FAIL = "\033[91m"
    RESET = "\x1b[0m"


# Logging configuration
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")


def resolve_ssh_config(host_alias, config_path="~/.ssh/config"):
    """
    Parse SSH configuration to resolve host details from ~/.ssh/config.

    Args:
        host_alias (str): The alias defined in SSH config (e.g., 'myvm').
        config_path (str): Path to the SSH config file.

    Returns:
        dict: Parsed configuration containing hostname, user, port, key, etc.
    """
    config_path = os.path.expanduser(config_path)
    ssh_config = paramiko.SSHConfig()
    with open(config_path) as f:
        ssh_config.parse(f)
    return ssh_config.lookup(host_alias)


def ssh_exec(client, command, alias=None):
    """
    Execute a command over an active SSH client connection.

    Args:
        client (paramiko.SSHClient): Active SSH connection.
        command (str): Command to run.
        alias (str, optional): Alias for display/logging.

    Returns:
        tuple: (stdout, stderr) as decoded strings.
    """
    _, stdout, stderr = client.exec_command(command)
    out, err = stdout.read().decode().strip(), stderr.read().decode().strip()

    if out:
        print(f"{bcolors.OKGREEN}[{alias or ''}: stdout] {out}{bcolors.RESET}")
    if err:
        print(f"{bcolors.FAIL}[{alias or ''}: stderr] {err}{bcolors.RESET}")

    return out, err


def ssh_exec_from_config(alias, command):
    """
    Connect via SSH using a config alias and run a command.

    Args:
        alias (str): The SSH alias to connect to (from ~/.ssh/config).
        command (str): The command to run.

    Returns:
        tuple: (stdout, stderr) as decoded strings.
    """
    info = resolve_ssh_config(alias)
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    client.connect(
        hostname=info["hostname"],
        username=info.get("user"),
        port=int(info.get("port", 22)),
        key_filename=info.get("identityfile", [None])[0],
    )

    try:
        return ssh_exec(client, command, alias)
    finally:
        client.close()


if __name__ == "__main__":
    import tkinter as tk
    from tkinter import simpledialog

    root = tk.Tk()
    root.withdraw()

    host_alias = "GPU-PC"
    remote_script = simpledialog.askstring(
        "Remote Script Path", "Enter path to remote script (e.g., ~/filter_bam.sh):"
    )
    input_dir = simpledialog.askstring(
        "Remote Input Dir", "Enter remote input directory path:"
    )
    output_dir = simpledialog.askstring(
        "Remote Output Dir", "Enter remote output directory path:"
    )

    # Build and execute command remotely
    command = f"ls '{input_dir}' '{output_dir}'"
    ssh_exec_from_config(host_alias, command)
