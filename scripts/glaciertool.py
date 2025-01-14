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


def check_metadata(bucket, key):
    """Fetch and print metadata for a single file."""
    try:
        metadata = s3_client.head_object(Bucket=bucket, Key=key)
        print(f"Metadata for {key}:\n{metadata}\n")
    except ClientError as e:
        print(f"Error retrieving metadata for {key}: {e}")


def check_thaw_status(bucket, key):
    """Check Glacier thaw status for a single file."""
    try:
        metadata = s3_client.head_object(Bucket=bucket, Key=key)
        storage_class = metadata.get("StorageClass", "STANDARD")
        restore_status = metadata.get("Restore", "Not archived")
        print(
            f"{key} - Storage Class: {storage_class}, Restore Status: {restore_status}"
        )
    except ClientError as e:
        print(f"Error checking thaw status for {key}: {e}")


def thaw_file(bucket, key, tier="Standard", days=7):
    """Restore a Glacier file to make it temporarily available."""
    try:
        s3_client.restore_object(
            Bucket=bucket,
            Key=key,
            RestoreRequest={"Days": days, "GlacierJobParameters": {"Tier": tier}},
        )
        print(f"Restore request initiated for {key} (Tier: {tier}, Days: {days})")
    except ClientError as e:
        print(f"Error restoring {key}: {e}")


def main():
    parser = argparse.ArgumentParser(description="S3 Glacier Management Tool")

    # Add argument for S3 URI first
    parser.add_argument("uri", help="S3 URI (e.g., s3://bucket/prefix/).")

    # Add argument for the action second (list, metadata, thaw-status, thaw)
    parser.add_argument(
        "action",
        choices=["list", "metadata", "thaw-status", "thaw"],
        help="Action to perform: list, metadata, thaw-status, or thaw.",
    )

    # Add optional arguments for tier and days
    parser.add_argument(
        "--tier", default="Standard", help="Thaw tier (Expedited/Standard/Bulk)."
    )
    parser.add_argument(
        "--days", type=int, default=7, help="Number of days to keep files available."
    )

    args = parser.parse_args()

    # Parse S3 URI
    bucket, prefix = parse_s3_uri(args.uri)

    if args.action == "list":
        files = list_files(bucket, prefix)
        print("Files found:")
        for file in files:
            print(file)
    else:
        files = list_files(bucket, prefix)
        if not files:
            print("No files found.")
            return

        print(f"Found {len(files)} file(s) under the prefix.")
        for file in files:
            if args.action == "metadata":
                check_metadata(bucket, file)
            elif args.action == "thaw-status":
                check_thaw_status(bucket, file)
            elif args.action == "thaw":
                thaw_file(bucket, file, args.tier, args.days)


if __name__ == "__main__":
    main()
