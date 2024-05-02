#!/bin/bash
# A simple script to manage plugins as git submodules (originally for vim plugins)

function vimfect() {
    cd "$(dirname "${BASH_SOURCE[0]}")"             || { echo "Failed to change directory to $script_dir"; exit 1; }
    git rev-parse --is-inside-work-tree &>/dev/null || { echo "Error: not inside a Git repository." && return 1; }
    [[ $# -eq 0 ]] && { echo "$0 is syncing plugins..."; git submodule sync && git submodule update --init --recursive; return; }
    [[ $# -eq 1 ]] && { echo "$0 added $1 submodule..."; git submodule add "$1" && git commit -m "add $1 as submodule"; return; }
}

vimfect "$@"

