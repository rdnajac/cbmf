#!/bin/bash
# Find .md5 files and check checksums concurrently without overloading the cores

find . -name "*.md5" -print0 | xargs -0 -n 1 -P "$(nproc)" md5sum -c
