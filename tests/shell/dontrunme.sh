#!/bin/bash
#
## Tests sourcing a script against running it.


SUCCESS()
{
  printf "\e[32mS\e[33mu\e[34mc\e[35cc\e[36me\e[37ms\e[0m!✨\n"
}

FAILURE()
{
  printf "\e[31mF\e[33ma\e[34mi\e[35ml\e[36me\e[37md\e[0m!💥\n"
}

main()
{
	printf "Hello, World!\n"
}

# Determine if the script is being sourced or executed.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
	__name__="__main__"
else
	__name__="__source__"
fi

# Code entry point. Only run if this script is being **run**, NOT sourced
# https://stackoverflow.com/a/70662116/4561887
if [ "$__name__" = "__main__" ]; then
	main "$@"
	printf "This script is being run, not sourced.\n"
	FAILURE
elif [ "$__name__" = "__source__" ]; then
	SUCCESS
else
	printf "Unknown state: %s\n" "$__name__"
	exit 1
fi
