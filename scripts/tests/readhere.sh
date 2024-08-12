#!/bin/bash

# Test read here document
# This script will call ./read.sh and feed it a heredoc
echo "This is a test"
./read.sh <<EOF
~/cbmf
EOF
echo "This is the end of the test"
