# from https://github.com/tmux-python/tmuxp/blob/master/src/tmuxp/util.py
import os
import pathlib
import shlex
import subprocess
import sys
import typing as t


def run_before_script(
    script_file: t.Union[str, pathlib.Path],
    cwd: t.Optional[pathlib.Path] = None,
) -> int:
    """Execute a shell script, wraps :meth:`subprocess.check_call()` in a try/catch."""
    try:
        proc = subprocess.Popen(
            shlex.split(str(script_file)),
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            cwd=cwd,
        )
        if proc.stdout is not None:
            for line in iter(proc.stdout.readline, b""):
                sys.stdout.write(console_to_str(line))
        proc.wait()

        if proc.returncode and proc.stderr is not None:
            stderr = proc.stderr.read()
            proc.stderr.close()
            stderr_strlist = console_to_str(stderr).split("\n")
            stderr_str = "\n".join(list(filter(None, stderr_strlist)))  # filter empty

            raise exc.BeforeLoadScriptError(
                proc.returncode,
                os.path.abspath(script_file),  # NOQA: PTH100
                stderr_str,
            )
    except OSError as e:
        if e.errno == 2:
            raise exc.BeforeLoadScriptNotExists(
                e,
                os.path.abspath(script_file),  # NOQA: PTH100
            ) from e
        raise
    return proc.returncode
