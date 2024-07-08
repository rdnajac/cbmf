#!/bin/bash
#
## Add a submodule to the repository.

set -euo pipefail
set -x

warn() {
	printf "\e[31m%s\e[0m\n" "$*"
}

bail() {
	warn "$*"
	exit 1
}

# Check if the user has provided a github url
[[ -n "$1" ]] || bail "Usage: $0 <github_url>"

if [[ "$1" =~ ^https://github.com ]]; then
	github_url="$1"
	
	# if there's a ? in the url, remove it and everything after it
	if [[ "$github_url" == *"?"* ]]; then
		github_url="${github_url%%\?*}"
	fi
else
	github_url="https://github.com/${1}.git"
fi

cd "${HOME}/cbmf/src" || exit 2

git submodule add --progress --force --depth 1 "$github_url"
