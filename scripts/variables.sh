#!/usr/bin/env sh

MMM_ROOT="$HOME/micromamba"
CURRENT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_ROOT=$(cd "$CURRENT_DIR/.." && pwd)
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
UTILS_DIR="$PROJECT_ROOT/utils"
TESTS_DIR="$PROJECT_ROOT/tests"

MAMBA_ROOT_PREFIX="$MMM_ROOT"
MAMBA_EXE="$MMM_ROOT/bin/micromamba"

ENVIRONMENTS_DIR="$PROJECT_ROOT/environments"

DEFAULT_ENV="cbmf"
MAMBA_INSTALL_URL="https://micro.mamba.pm/api/micromamba"
