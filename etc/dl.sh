#!/bin/sh
#===============================================================================
## @brief Download a binary from a URL and make it executable
## @param $1 name of the binary
## @param $2 URL of the binary
## @return exit code 0 if successful
#===============================================================================
download_binary() {
	if [ "$#" -ne 2 ]; then
		err_exit "Usage: download_binary <name> <url>" EINVAL
	fi

	name="$1"
	url="$2"
	mkdir -p "$BIN_FOLDER"


	if hash curl > /dev/null 2>&1; then
		curl "$url" -o "${BIN_FOLDER}/${name}" -fsSL --compressed "${CURL_OPTS:-}"
	elif hash wget > /dev/null 2>&1; then
		wget "$url" -O "${BIN_FOLDER}/${name}" -q --show-progress
	else
		err_exit "Neither curl nor wget were found" ENOCMD
	fi
	chmod +x "${BIN_FOLDER}/${name}"
}

