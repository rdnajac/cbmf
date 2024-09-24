#!/usr/bin/env sh
set -eux

CURRENT_DIR=$(cd $(dirname "$0") && pwd)
UTILS_DIR="${CURRENT_DIR}/../utils"

. "$UTILS_DIR/utils.sh"
. "$UTILS_DIR/install_micromamba.sh"

check_prereqs
detect_platform_architecture
detect_shell
install_micromamba $HOME


