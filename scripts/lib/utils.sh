#!/bin/bash
#
# Source this file to configure the script environment
# https://gabrielstaples.com/bash-libraries
# https://superuser.com/questions/46139/what-does-source-do
# https://askubuntu.com/questions/182012/is-there-a-difference-between-and-source-in-bash-after-all

# Define ANSI color codes for colored output
export RED='\033[0;31m'
export GRN='\033[0;32m'
export YEL='\033[0;33m'
export BLU='\033[0;34m'
export MAG='\033[0;35m'
export CYN='\033[0;36m'
export WHT='\033[0;37m'
export END='\033[0m'

# Helper functions to print colored messages
warn()
{
	printf "${RED}%s${RESET}\n" "$1" >&2
}

info()
{
	printf "${BLU}%s${RESET}\n" "$1" >&2
}

# Quit with an error message and optional exit status
bail()
{
	warn "$1"
	exit "${2:-1}"
}

assert()
{
	local condition="$1"
	shift
	[[ "$condition" ]] || bail "$@"
}

# Error Codes
readonly E_COMMAND_NOT_EXECUTABLE=126
readonly E_COMMAND_NOT_FOUND=127

# Check if a command is available
ensure_installed()
{
	# [[ -e "$1" ]] || bail "Error: $1 not found" "$E_COMMAND_NOT_FOUND"
	# check another way
	[[ -x "$(command -v "$1")" ]] || bail "Error: $1 not found" "$E_COMMAND_NOT_FOUND"

	[[ -x "$1" ]] || bail "Error: $1 is not executable" "$E_COMMAND_NOT_EXECUTABLE"
}

# Test
_test()
{
	warn "Running tests... from ${BASH_SOURCE[0]}"
	ensure_installed "ls" && info "ls is installed and executable"

	assert [[ -f "$0" ]] "This test should pass"

	info "All tests passed"
}

main()
{
	_test
}

# Determine if the script is being sourced or executed (run).
# See:
if [ """${BASH_SOURCE[0]}" = """$0" ]; then
	# This script is being run.
	__name__="__main__"
else
	# This script is being sourced.
	__name__="__source__"
fi

# Code entry point. Only run if this script is being **run**, NOT sourced
# https://stackoverflow.com/a/70662116/4561887
if [ "$__name__" = "__main__" ]; then
	main "$@"
fi
