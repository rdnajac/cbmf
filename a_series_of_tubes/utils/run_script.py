import subprocess
import pathlib
import shlex
import os
from typing import Union, Optional
from .logger import logger


def run_script(
    script_file: Union[str, pathlib.Path],
    cwd: Optional[pathlib.Path] = None,
) -> int:
    """
    Execute a shell script and handle errors.

    Args:
        script_file: Path to the script file to be executed.
        cwd: Directory to change to before executing the script.

    Returns:
        The return code of the script execution.

    Raises:
        subprocess.CalledProcessError: If the script exits with a non-zero status.
        FileNotFoundError: If the script file does not exist.
    """
    # Expand the user's home directory if the path starts with ~
    script_file_str = os.path.expanduser(str(script_file))

    try:
        logger.info(f"Executing script: {script_file_str}")
        result = subprocess.run(
            shlex.split(script_file_str),
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        print(result.stdout)
        return result.returncode
    except FileNotFoundError:
        print(f"Script file {script_file_str} does not exist")
        logger.error(f"Script file {script_file_str} does not exist")
        raise
    except subprocess.CalledProcessError as e:
        print(f"Script execution failed with return code {e.returncode}")
        print(f"Error output:\n{e.stderr}")
        logger.error(f"Script execution failed with return code {e.returncode}")
