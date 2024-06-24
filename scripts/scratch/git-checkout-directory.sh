#!/bin/bash

# Sources:
# https://web.archive.org/web/20210309072829/http://jasonkarns.com/blog/2011/11/15/subdirectory-checkouts-with-git-sparse-checkout/
# https://stackoverflow.com/questions/2416815/how-to-git-pull-all-but-one-folder/17075665#17075665#
# https://stackoverflow.com/questions/2425059/how-to-pull-specific-directory-with-git

set -euo pipefail

git_checkout_directory() {
  local uri="$1"

  # Strip out '/tree/master' if present and remove trailing slash
  local sanitized_url="${uri%%/tree/master*}"
  sanitized_url="${sanitized_url%/}"

  # Validate the sanitized URL
  if ! [[ "$sanitized_url" =~ ^https:\/\/github\.com\/[^\/]+\/[^\/]+\/?$ ]]; then
    echo "Invalid URL: $sanitized_url" >&2
    return 1
  fi

  # Get substring after 'tree/master/' and strip trailing slash if exists
  local directory="${uri##*tree/master/}"
  directory="${directory%/}"

  # Extract the final part of the path to use as the target directory
  local target_dir="${directory##*/}"
  mkdir -p "$target_dir"
  cd "$target_dir" || {
    echo "Failed to create or change to directory $target_dir." >&2
    return 1
  }

  git init
  git remote add -f origin "$sanitized_url"
  git config core.sparsecheckout true
  echo "$directory/" >> .git/info/sparse-checkout
  git pull origin master

  # Safely move contents from the nested directory to the target directory
  find "$directory" -mindepth 1 -maxdepth 1 -exec mv {} . \;

  # Optionally remove the empty directory and the .git directory
  rmdir "$directory"
  rm -rf .git

  echo "Done. Repository folder copied to $target_dir."
}

main() {
  if [[ -z "$1" ]]; then
    echo "Usage: $0 <repository-url>"
    exit 1
  fi

  git_checkout_directory "$1"
}

main "$@"
