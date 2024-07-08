#!/bin/bash

## Add a git submodule to the repository
set -euo pipefail
set -x

# Check if the user has provided a github url
if [ "$1" = "" ]; then
	echo "Please provide a github url"
	exit 1
fi

github_url=$1

# if the githuburl doesnt begin with https://github.com, then add it
if [[ ! "$github_url" =~ ^https://github.com ]]; then
	github_url="https://github.com/${github_url}.git"
fi

PLUGIN_DIR="${HOME}/.vim/pack/vimfect/opt"
mkdir -p "$PLUGIN_DIR"
cd "$PLUGIN_DIR" || exit 1

git submodule add --depth=1 --force "$github_url"
