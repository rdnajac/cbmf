#!/usr/bin/env sh

set -eu

THIS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "$0")" && pwd)}"

. "$THIS_DIR/utils.sh"
. "$THIS_DIR/variables.sh"

# Check for required commands
for cmd in tar bzip2 curl; do
	command -v "$cmd" > /dev/null 2>&1 || err_exit "$cmd is required but not found"
done

if [ ! -f "$MAMBA_EXE" ]; then
	# Determine platform and architecture
	PLATFORM_ARCH=$(detect_platform_architecture)
	RELEASE_URL="https://micro.mamba.pm/api/micromamba/$PLATFORM_ARCH/latest"

	# Download and extract micromamba
	mkdir "$HOME"/micromamba && cd "$HOME"/micromamba
	curl -Ls "$RELEASE_URL" | tar -xvj bin/micromamba
	chmod +x "$MAMBA_EXE" || true # Ignore errors
	printf "Downloaded micromamba to %s\n" "$MAMBA_EXE"
	"$MAMBA_EXE" shell init -s "${SHELL##*/}" -p "$MAMBA_ROOT_PREFIX"
	printf "Initialized micromamba shell settings. Restarting shell...\n"
	exec "$SHELL"
fi

# initialize_micromamba_shell "$MAMBA_ROOT_PREFIX" "$MAMBA_EXE"
# write_condarc
eval "$("$MAMBA_EXE" shell hook -s "$DETECTED_SHELL")" || err_exit "Failed to set up micromamba shell hook"
SUCCESS
cbmf_default_env
