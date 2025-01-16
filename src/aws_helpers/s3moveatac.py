#!/usr/bin/env python3

import boto3
import subprocess
from urllib.parse import urlparse

# Initialize S3 client
s3_client = boto3.client("s3")


def parse_s3_uri(uri):
    """Parse an S3 URI into bucket and prefix."""
    parsed = urlparse(uri)
    if parsed.scheme != "s3":
        raise ValueError(f"Invalid S3 URI: {uri}")
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")
    return bucket, prefix


def list_subfolders(bucket, prefix):
    """List all subfolders within a given prefix."""
    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        operation_params = {"Bucket": bucket, "Prefix": prefix, "Delimiter": "/"}
        subfolders = []
        for page in paginator.paginate(**operation_params):
            if "CommonPrefixes" in page:
                for common_prefix in page["CommonPrefixes"]:
                    subfolders.append(common_prefix["Prefix"])
        return subfolders
    except Exception as e:
        print(f"Error listing subfolders: {e}")
        return []


def main():
    source_uri = "s3://lab-aaf-ngs-data-archive/ATACseq/20180711_Dx-Re_ATACseq_FG/"
    dest_uri = "s3://lab-aaf-ngs-data-archive/ATACseq/20180711_ATACseq_Dx-Re_FG/"

    # Parse the source S3 URI
    source_bucket, source_prefix = parse_s3_uri(source_uri)

    # List all subfolders
    subfolders = list_subfolders(source_bucket, source_prefix)

    if not subfolders:
        print(f"No subfolders found in {source_uri}")
        return

    # Call s3move.py for each subfolder
    for subfolder in subfolders:
        full_source = f"s3://{source_bucket}/{subfolder}"
        print(f"Processing subfolder: {full_source}")
        subprocess.run(["./s3move.py", full_source, dest_uri], check=True)


if __name__ == "__main__":
    main()
