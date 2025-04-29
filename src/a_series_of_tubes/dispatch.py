#!/usr/bin/env python3
# Dispatch shell commands with timestamped logging

import subprocess
import sys
import time
import shlex


def log(msg):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}", file=sys.stderr)


def run(cmd):
    log(f"Running: {cmd}")
    start = time.time()
    try:
        subprocess.run(shlex.split(cmd), check=True)
    except subprocess.CalledProcessError as e:
        log(f"Command failed with exit code {e.returncode}")
        sys.exit(e.returncode)
    duration = time.time() - start
    log(f"Finished in {duration:.2f}s")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} '<command>'", file=sys.stderr)
        sys.exit(1)
    cmd = " ".join(sys.argv[1:])
    run(cmd)
