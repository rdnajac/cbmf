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
    curl -Ls "$RELEASE_URL" | tar -xvj bin/micromamba
    chmod +x "$MAMBA_EXE"
    printf "Downloaded micromamba to %s\n" "$MAMBA_EXE"
    "$MAMBA_EXE" shell init -s "${SHELL##*/}" -p "$MAMBA_ROOT_PREFIX"
    printf "Initialized micromamba shell settings. Restarting shell...\n"
    exec "$SHELL"
fi

# Initialize micromamba shell
initialize_micromamba_shell "$MAMBA_ROOT_PREFIX" "$MAMBA_EXE"

# Write .condarc file
write_condarc

# Create default environment
cbmf_default_env
