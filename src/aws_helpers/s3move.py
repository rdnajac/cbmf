#!/usr/bin/env python3

import boto3
import argparse
from urllib.parse import urlparse
from botocore.exceptions import ClientError

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


def move_file(source_bucket, source_key, dest_bucket, dest_key):
    """
    Move a file from one S3 location to another.
    
    Args:
        source_bucket (str): The bucket containing the source file.
        source_key (str): The key of the source file.
        dest_bucket (str): The destination bucket.
        dest_key (str): The key for the destination file.
    """
    try:
        # Copy the file to the new location
        s3_client.copy_object(
            CopySource={"Bucket": source_bucket, "Key": source_key},
            Bucket=dest_bucket,
            Key=dest_key,
        )
        print(f"Copied {source_key} to {dest_bucket}/{dest_key}")

        # Delete the original file
        s3_client.delete_object(Bucket=source_bucket, Key=source_key)
        print(f"Deleted {source_key} from {source_bucket}")
    except ClientError as e:
        print(f"Error moving {source_key}: {e}")


def list_files(bucket, prefix):
    """List all files in the specified bucket and prefix using boto3."""
    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        operation_params = {"Bucket": bucket, "Prefix": prefix}
        file_list = []
        for page in paginator.paginate(**operation_params):
            if "Contents" in page:
                for obj in page["Contents"]:
                    file_list.append(obj["Key"])
        return file_list
    except ClientError as e:
        print(f"Error listing files: {e}")
        return []


def main():
    parser = argparse.ArgumentParser(description="Move files between S3 locations")
    parser.add_argument("source", help="Source S3 URI (e.g., s3://bucket/prefix/).")
    parser.add_argument("destination", help="Destination S3 URI (e.g., s3://bucket/prefix/).")

    args = parser.parse_args()

    # Parse S3 URIs
    source_bucket, source_prefix = parse_s3_uri(args.source)
    dest_bucket, dest_prefix = parse_s3_uri(args.destination)

    # List all files in the source
    files = list_files(source_bucket, source_prefix)
    if not files:
        print(f"No files found in {args.source}")
        return

    print(f"Found {len(files)} file(s) in {args.source}. Starting move process...")

    for file in files:
        # Compute destination key
        dest_key = f"{dest_prefix}{file[len(source_prefix):]}"
        move_file(source_bucket, file, dest_bucket, dest_key)


if __name__ == "__main__":
    main()
