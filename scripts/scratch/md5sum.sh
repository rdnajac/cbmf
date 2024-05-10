# Generate checksums for all files in the current directory
generate_checksums()
{
    local extension=${1:-"fastq.gz"}
    local checksum_file="md5sums.txt"
    echo "Generating checksums for *.$extension files in the current directory..."
    for file in *.$extension; do
        echo "$(md5sum "$file" | cut -d ' ' -f 1)\t$file" >> $checksum_file
    done
    echo "Checksums stored in $checksum_file"
}

# Verify checksums for all files in the current directory
verify_checksums()
{
    local md5sums_txt="${1:-md5sums.txt}"
    echo "Verifying checksums..."
    while read -r expected_checksum file; do
        local status="[ OK ]"
        [[ "$(md5sum "$file" | cut -d ' ' -f 1)" != "$expected_checksum" ]] && status="[FAIL]"
        echo "${status} ${file}"
    done < $md5sums_txt
}

