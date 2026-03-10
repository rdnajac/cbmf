#!/usr/bin/env python3
import sys
from remote_runner import execute_remote_command
from logger import logger

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: ./main.py <ssh-alias> <command>")
        sys.exit(1)
    alias, cmd = sys.argv[1], " ".join(sys.argv[2:])
    try:
        execute_remote_command(alias, cmd)
    except Exception as e:
        logger.exception("Execution failed")
        sys.exit(1)
