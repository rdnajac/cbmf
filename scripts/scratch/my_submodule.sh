#!/bin/bash

set -euxo pipefail

mv_submodule() {
    local old_path=$1
    local new_path=$2
    local submodule_url

    if [[ -z "$old_path" || -z "$new_path" ]]; then
        echo "Usage: mv_submodule <old_submodule_path> <new_submodule_path>"
        return 1
    fi

    if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
        echo "Error: Not inside a Git repository"
        return 1
    fi

    submodule_url=$(git config -f .gitmodules --get "submodule.${old_path}.url")
    if [[ -z "$submodule_url" ]]; then
        echo "Error: Submodule URL not found for $old_path"
        return 1
    fi

    git submodule deinit -f "$old_path"
    rm -rf ".git/modules/$old_path"
    git rm -f "$old_path"

    git submodule add "$submodule_url" "$new_path"

    git mv .gitmodules .gitmodules.bak
    sed "s|$old_path|$new_path|g" .gitmodules.bak > .gitmodules
    rm .gitmodules.bak
    git add .gitmodules

    git commit -m "Moved submodule from $old_path to $new_path"
}

# Example usage
# mv_submodule libs/my_submodule src/my_submodule
