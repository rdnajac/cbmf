import paramiko

def resolve_ssh_config(host_alias, config_path="~/.ssh/config"):
    config_path = os.path.expanduser(config_path)
    ssh_config = paramiko.SSHConfig()
    with open(config_path) as f:
        ssh_config.parse(f)
    return ssh_config.lookup(host_alias)

def connect_ssh(alias):
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
    return client

def ssh_exec(client, command):
    _, stdout, stderr = client.exec_command(command)
    return stdout.read().decode().strip(), stderr.read().decode().strip()
