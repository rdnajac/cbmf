#!/bin/bash
# ==============================================================================
# Usage: ./scriptname.sh -g -v [checksum_file] -x [extension]
# -g: Generate a new checksum file from files of the specified extension in the current directory.
# -v: Verify checksums using a provided checksum file.
# -x: Specify the file extension to look for (default is 'fastq.gz').
# ==============================================================================

extension="fastq.gz"  # Default file extension

function generate_checksums {
    echo "Generating checksums for *.$extension files in the current directory..."
    rm -f "$checksum_file"  # Remove existing checksum file if present
    for file in *.$extension; do
        echo "$(md5sum "$file" | cut -d ' ' -f 1) $file" >> "$checksum_file"
    done
    echo "Checksums stored in $checksum_file"
}

function verify_checksums {
    echo "Verifying checksums using $checksum_file..."
    while read -r expected_checksum file; do
        [[ ! -f "$file" ]] && echo "File not found: $file" && continue
        [[ "$(md5sum "$file" | cut -d ' ' -f 1)" == "$expected_checksum" ]] && echo "$file: Checksums match" || echo "$file: Checksums do not match"
    done < "$checksum_file"
}

while getopts ":g:v:x:" opt; do
    case "$opt" in
        g)
            checksum_file="$OPTARG"
            action="generate"
            ;;
        v)
            checksum_file="$OPTARG"
            action="verify"
            ;;
        x)
            extension="$OPTARG"
            ;;
        *)
            echo "Usage: $0 [-g|-v checksum_file] [-x extension]"
            exit 1
            ;;
    esac
done

[[ -z "$checksum_file" ]] && echo "Checksum file not specified." && exit 1
[[ "$action" == "generate" ]] && generate_checksums && exit 0
[[ "$action" == "verify" ]] && verify_checksums && exit 0
echo "Invalid or no action specified."
echo "Usage: $0 [-g|-v checksum_file] [-x extension]"
exit 1

