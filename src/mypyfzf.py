#!/usr/bin/env python3

import os
import subprocess
from shutil import which

FZF_URL = "https://github.com/junegunn/fzf"


class FzfPicker:
    def __init__(self, executable_path=None):
        if executable_path:
            self.executable_path = executable_path
        elif not which("fzf"):
            raise SystemError(f"'fzf' is not found in PATH. Install it from {FZF_URL}")
        else:
            self.executable_path = "fzf"

    def _run_fzf(self, input_lines, header=None, fzf_options=""):
        """Run fzf with given input and options, return selected item or None."""
        cmd = [self.executable_path]
        if header:
            cmd.extend(["--header", header])
        if fzf_options:
            cmd.extend(fzf_options.split())

        try:
            result = subprocess.run(
                cmd,
                input="\n".join(str(line) for line in input_lines),
                capture_output=True,
                text=True,
                check=True,
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            return None

    def pick(self, items, header=None, fzf_options=""):
        """Pick a single item from a list."""
        return self._run_fzf(items, header, fzf_options)

    def pick_dir(self, start_path=".", header="Select a directory", fzf_options=""):
        """Pick a directory from the filesystem using fd (if available) or find."""
        fd_cmd = "fd" if which("fd") else ("fdfind" if which("fdfind") else None)
        try:
            if fd_cmd:
                list_dirs_cmd = [
                    fd_cmd,
                    "--type",
                    "d",
                    "--hidden",
                    "--follow",
                    "--exclude",
                    ".git",
                    ".",
                    start_path,
                ]
            else:
                list_dirs_cmd = ["find", start_path, "-type", "d"]

            result = subprocess.run(
                list_dirs_cmd, capture_output=True, text=True, check=True
            )
            dirs = result.stdout.strip().splitlines()
            return self._run_fzf(dirs, header, fzf_options)
        except subprocess.CalledProcessError:
            return None


if __name__ == "__main__":
    fzf = FzfPicker()
    home_dir = os.path.expanduser("~")  # Expand ~ into /home/username
    selection = fzf.pick_dir(start_path=home_dir)
    print("You selected:", selection)
