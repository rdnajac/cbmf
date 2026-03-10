import subprocess
import pathlib
import shlex
import os
from logger import logger

def run_script(script_file, cwd=None):
    path = os.path.expanduser(str(script_file))
    if not os.path.exists(path):
        logger.error(f"Script file {path} does not exist")
        raise FileNotFoundError(f"{path} not found")

    logger.info(f"Executing script: {path}")
    result = subprocess.run(
        shlex.split(path),
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    logger.info(f"Script output:\n{result.stdout}")
    if result.stderr:
        logger.warning(f"Script stderr:\n{result.stderr}")
    if result.returncode != 0:
        logger.error(f"Script failed with return code {result.returncode}")
        raise subprocess.CalledProcessError(result.returncode, path)
    return result.returncode