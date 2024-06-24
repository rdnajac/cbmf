#!/bin/bash
#
## Quick and dirty script to clone a git repo, install it and add it to the PATH

function clone_cd_make_export
{
	git clone "$1"
	cd "$(basename "$1" .git)" || exit 1
	make -j "$(nproc)"
	export PATH=$PATH:$PWD
}

clone_cd_make_export "$1"
