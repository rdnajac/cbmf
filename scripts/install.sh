#!/bin/bash
# cross-platform installation script

set -eu # Safety first!
# This is the same invoking both `set -o errexit` and `set -o nounset`.
# It causes the script to exit immediately if a command exits with a 
# non-zero status, and treats unset variables error.


# Catch the error in case a pipeline element fails (returns non-zero)
# Note that this does not catch every error, but a script that works
# with `set -o pipefail` will always work without it, while the opposite is not true.
# Note that `set -o pipefail` is not POSIX-compliant, so we can add a check for POSIX mode.
if [ "${POSIXLY_CORRECT:-}" != "" ]; then
	printf "POSIXLY_CORRECT is set; some features may be disabled\n"
else
	set -o pipefail
fi

# Uncomment to print each command before it is executed. Useful for debugging.
# This is the same as `set -o xtrace` or invoking the script with `bash -x`.
set -x

# ==============================================================================
# Constants
# ==============================================================================
# The 0th positional parameter is the name of the script itself
# the cd && pwd ensure that the script works even if it is called 
# from a relative directory, symlink, or over ssh
CBMF_DIR=$(cd "$(dirname "$0")" && pwd)
export BIN_DIR="${CBMF_DIR}/bin"
export DEV_DIR="${CBMF_DIR}/dev"
export GENOMES_DIR="${CBMF_DIR}/genomes"
export SCRIPTS_DIR="${CBMF_DIR}/scripts"

mkdir -p "$BIN_DIR" "$GENOMES_DIR" "$SCRIPTS_DIR"

export MAMBA_ROOT_PREFIX="${HOME}/cbmf/micromamba"
export MAMBA_EXE="${HOME}/cbmf/bin/micromamba"

# This is only env variable relevant if you are using R with reticulate,
# an R package that allows you to run Python code in R.
export RETICULATE_MINICONDA_ENABLED=FALSE
# Entering the following commands in R should NOT prompt you to install miniconda:
# > library(reticulate)
# > os <- import("os")

#===============================================================================
# Utility functions
#===============================================================================
err() { printf "Error: %s\n" "$1" >&2; }

err_exit() { err "$1"; exit "${2:-1}"; }

#===============================================================================
# Check for required commands
# ==============================================================================
for cmd in tar bzip2 curl; do
	if ! command -v "$cmd" > /dev/null 2>&1; then
		err_exit "$cmd is required but not found"
	fi
done

# ==============================================================================
# Determine platform and architecture
# ==============================================================================
case "$(uname)-$(uname -m)" in
	Linux-x86_64)  PLATFORM_ARCH="linux-64" ;;
	Linux-aarch64) PLATFORM_ARCH="linux-aarch64" ;;
	Linux-ppc64le) PLATFORM_ARCH="linux-ppc64le" ;;
	Darwin-x86_64) PLATFORM_ARCH="osx-64" ;;
	Darwin-arm64)  PLATFORM_ARCH="osx-arm64" ;;
	*) err_exit "Unsupported platform-architecture combination" ;;
esac

# ==============================================================================
# Download and extract micromamba 
#
# The following magic URL always returns the latest available version of micromamba,
# and the bin/micromamba part is automatically extracted using tar.
# ==============================================================================
RELEASE_URL="https://micro.mamba.pm/api/micromamba/$PLATFORM_ARCH/latest"

# Create directory, download and extract micromamba
# If we already have a micromamba executable, we skip this step
curl -Ls "$RELEASE_URL" | tar -xvj bin/micromamba
chmod +x "$MAMBA_EXE"
printf "Downloaded micromamba to %s\n" "$MAMBA_EXE"
./bin/micromamba shell init -s "${SHELL##*/}" -p "$MAMBA_ROOT_PREFIX"
printf "Initialized micromamba shell settings. Restarting shell...\n"
exec "$SHELL"
# $MAMBA_EXE shell init --shell bash --root-prefix "$MAMBA_ROOT_PREFIX"


# ==============================================================================
# micromamba configuration
#
# These commands create a ~/.condarc file with the above settings.
# The ~/.condarc file is read by micromamba and conda, and can be edited manually.
# However, CLI options take precedence over settings in the ~/.condarc file,
# so a more reliable way to set these options is to use wrapper scripts to
# manually set the options before calling micromamba.
#
# Rather, just use the .yaml file to specify the channels and channel_priority.
# ==============================================================================
# micromamba config append channels bioconda
# micromamba config append channels conda-forge
# micromamba config append channels nodefaults
# micromamba config set channel_priority strict

# ==============================================================================
# Install the default environment
# ==============================================================================
# micromamba create -n cbmf -f "${DEV_DIR}/cbmf.yml"
# micromamba create -n analysis -f ./dev/analysis.yml
# micromamba activate cbmf

# RStudio
# TODO: if [ -n "$RSTUDIO" ] or if [ -n "$RSTUDIO_CONSOLE" ]
# note: this is how R sets the prompt: declare a variable, echo it, then unset it
# PROMPT_COMMAND=(_RS_PWD=$(dirs +0); echo -ne "\033]0;${_RS_PWD}\007"; unset _RS_PWD)

exit 0