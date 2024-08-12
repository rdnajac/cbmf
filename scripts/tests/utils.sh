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
export RESET='\033[0m'

# # Check if a command is available
# ensure_installed()
# {
# 	# [[ -e "$1" ]] || bail "Error: $1 not found" "$E_COMMAND_NOT_FOUND"
# 	# check another way
# 	[[ -x "$(command -v "$1")" ]] || bail "Error: $1 not found" "$E_COMMAND_NOT_FOUND"

# 	[[ -x "$1" ]] || bail "Error: $1 is not executable" "$E_COMMAND_NOT_EXECUTABLE"
# }

# # Test
# _test()
# {
# 	warn "Running tests... from ${BASH_SOURCE[0]}"
# 	ensure_installed "ls" && info "ls is installed and executable"

# 	assert [[ -f "$0" ]] "This test should pass"

# 	info "All tests passed"
# }

main()
{
	printf "Hello, World!\n"
}

# Determine if the script is being sourced or executed.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
	# This script is being run.
	__name__="__main__"
else
	# This script is being sourced.
	__name__="__source__"
	info "Sourced ${BASH_SOURCE[0]}"
fi

# Code entry point. Only run if this script is being **run**, NOT sourced
# https://stackoverflow.com/a/70662116/4561887
if [ "$__name__" = "__main__" ]; then
	main "$@"
fi
