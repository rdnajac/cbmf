import os
import posixpath
import ntpath
from ssh_utils import connect_ssh, ssh_exec
from logger import logger

def execute_remote_script(alias, local_path, remote_dir="~/remote_scripts"):
    local_path = os.path.expanduser(local_path)
    if not os.path.isfile(local_path):
        raise FileNotFoundError(local_path)

    client = connect_ssh(alias)
    sftp = client.open_sftp()

    remote_dir = ssh_exec(client, f"mkdir -p {remote_dir} && realpath {remote_dir}")[0]
    remote_name = ntpath.basename(local_path)
    remote_path = posixpath.join(remote_dir, remote_name)

    sftp.put(local_path, remote_path)
    sftp.chmod(remote_path, 0o755)
    sftp.close()

    logger.info(f"Transferred script to {alias}:{remote_path}")
    out, err = ssh_exec(client, f"{remote_path}")
    if out:
        logger.info(f"[{alias}] stdout:\n{out}")
    if err:
        logger.error(f"[{alias}] stderr:\n{err}")
    client.close()
    return out, err