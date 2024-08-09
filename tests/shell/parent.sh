#!/bin/sh
DETECTED_PARENT=$(ps -o comm "$PPID" | tail -1)
# Remove leading dash from login shells using the `#` syntax
# to remove the characters from the beginning of the string
DETECTED_PARENT=${DETECTED_PARENT#-}

printf "Hello, %s user!\n" "$DETECTED_PARENT"

case "$DETECTED_PARENT" in
	bash | fish | xonsh | zsh) shell=$DETECTED_PARENT ;;
	*) shell=${SHELL##*/} ;;
esac
# Fallback to $SHELL variable, which contains the path to the current
# shell. Trim the leading path from the variable to get the shell name.
# We can also set shell=$(basename "$SHELL")

