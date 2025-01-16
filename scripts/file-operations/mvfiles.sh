#!/bin/bash

# Define source and destination S3 URIs
SOURCE="s3://lab-aaf-ngs-data-archive/ATACseq/20180711_Dx-Re_ATACseq_FG/"
DESTINATION="s3://lab-aaf-ngs-data-archive/ATACseq/20180711_ATACseq_Dx-Re_FG/"

# List all folders in the source directory
folders=$(aws s3 ls "$SOURCE" --recursive | awk '{print $4}' | grep '/$')

# Loop through each folder
for folder in $folders; do
    echo "Processing folder: $folder"

    # Construct the full S3 paths
    full_source="${SOURCE}${folder}"
    full_destination="${DESTINATION}${folder}"

    # Sync all .gz files from this folder to the new destination
    aws s3 sync "$full_source" "$full_destination" --exclude "*" --include "*.gz" --dryrun

    echo "Synced .gz files from $full_source to $full_destination"
done

echo "All folders processed!"
