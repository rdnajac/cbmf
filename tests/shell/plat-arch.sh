#!/bin/sh

# uname returns (operating) system information
# We only support linux-64 (for now), but micromamba supports:
# linux-aarch64 | linux-ppc64le | linux-64 | osx-arm64 | osx-64)
# MICROMAMBA_URL="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-${PLATFORM}-${ARCHITECTURE}"

case "$(uname)" in
	Linux) PLATFORM="linux" ;;
	Darwin) PLATFORM="osx" ;;
esac

# uname -m returns the machine hardware name
ARCHITECTURE="$(uname -m)"
case "$ARCHITECTURE" in
	aarch64 | ppc64le | arm64) ;; # pass
	*) ARCHITECTURE="64" ;;
esac

case "$PLATFORM-$ARCHITECTURE" in
	linux-64 | osx-64) ;; # pass
	*) err_exit "Unsupported platform: $PLATFORM-$ARCHITECTURE" ;;
esac
