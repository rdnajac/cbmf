#!/bin/bash
#
## Main entry point for cbmf workflows. 
## Get options and call the appropriate script.

set -e  # Exit on error
set -u  # Exit on using unset variable
set -x  # Print commands

source lib/utils.sh

# Get options
while getopts "h" opt; do
  case $opt in
    h)
      echo "Usage: cbmf.sh [OPTIONS] COMMAND [ARGS]"
      echo "Options:"
      echo "  -h    Show this message and exit"
      echo "Commands:"
      echo "  run   Run the cbmf workflow"
      echo "  plot  Plot the results of the cbmf workflow"
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done
