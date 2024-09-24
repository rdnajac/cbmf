#!/bin/bash

# Function to print failure status
FAILURE() {
    printf "\e[31mF\e[33ma\e[34mi\e[35ml\e[36me\e[37md\e[0m!💥\n"
}

# Function to display current directory and script details
show_dir_info() {
    DIR_NO_CD=$(dirname "$0")
    DIR=$(cd "$(dirname "$0")" && pwd)

    echo "DIR_NO_CD: $DIR_NO_CD"
    echo "DIR: $DIR"
    echo "PWD: $PWD"
    echo "The script you are running has:"
    echo "basename: [$(basename "$0")]"
    echo "dirname : [$(dirname "$0")]"
    echo "pwd     : [$(pwd)]"
}

# Function to determine if the script is being sourced or executed
detect_script_execution() {
    if [ "${BASH_SOURCE[0]}" = "$0" ]; then
        __name__="__main__"
    else
        __name__="__source__"
    fi

    if [ "$__name__" = "__main__" ]; then
        main "$@"
        printf "This script is being run, not sourced.\n"
    elif [ "$__name__" = "__source__" ]; then
        FAILURE
    else
        printf "Unknown state: %s\n" "$__name__"
        exit 1
    fi
}

# Main function for demonstration
main() {
    printf "Hello, World!\n"
}

# Function to detect the current shell
detect_shell() {
    DETECTED_PARENT=$(ps -o comm "$PPID" | tail -1)
    DETECTED_PARENT=${DETECTED_PARENT#-}

    printf "Hello, %s user!\n" "$DETECTED_PARENT"

    case "$DETECTED_PARENT" in
        bash | fish | xonsh | zsh) shell=$DETECTED_PARENT ;;
        *) shell=${SHELL##*/} ;;
    esac

    echo "Current shell: $shell"
}

# Function to detect platform and architecture
detect_platform_architecture() {
    case "$(uname)" in
        Linux) PLATFORM="linux" ;;
        Darwin) PLATFORM="osx" ;;
        *) err_exit "Unsupported OS" ;;
    esac

    ARCHITECTURE="$(uname -m)"
    case "$ARCHITECTURE" in
        aarch64 | ppc64le | arm64) ;; # pass
        *) ARCHITECTURE="64" ;;
    esac

    case "$PLATFORM-$ARCHITECTURE" in
        linux-64 | osx-64) ;; # pass
        *) err_exit "Unsupported platform: $PLATFORM-$ARCHITECTURE" ;;
    esac
}

# Function to download and run an installation script with a heredoc
install_with_heredoc() {
    BIN_DIR=~/bin

    bash <(curl -L micro.mamba.pm/install.sh) << EOF
${BIN_DIR}
Y
n
EOF

    if [ -f "${BIN_DIR}/micromamba" ]; then
        echo "micromamba has been successfully installed in ${BIN_DIR}"
    else
        echo "Installation of micromamba may have failed. Please check ${BIN_DIR}"
    fi
}


