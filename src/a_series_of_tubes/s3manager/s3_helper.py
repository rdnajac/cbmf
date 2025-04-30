import boto3
from urllib.parse import urlparse
from botocore.exceptions import ClientError


class S3Helper:
    def __init__(self):
        self.client = boto3.client("s3")

    def parse_uri(self, uri):
        parsed = urlparse(uri)
        if parsed.scheme != "s3":
            raise ValueError(f"Invalid S3 URI: {uri}")
        return parsed.netloc, parsed.path.lstrip("/")

    def list_files(self, bucket, prefix):
        try:
            paginator = self.client.get_paginator("list_objects_v2")
            params = {"Bucket": bucket, "Prefix": prefix}
            return [
                obj["Key"]
                for page in paginator.paginate(**params)
                if "Contents" in page
                for obj in page["Contents"]
            ]
        except ClientError as e:
            print(f"Error listing files: {e}")
            return []

    def get_metadata(self, bucket, key):
        try:
            return self.client.head_object(Bucket=bucket, Key=key)
        except ClientError as e:
            print(f"Error retrieving metadata for {key}: {e}")

    def get_thaw_status(self, bucket, key):
        try:
            metadata = self.client.head_object(Bucket=bucket, Key=key)
            return metadata.get("StorageClass", "STANDARD"), metadata.get(
                "Restore", "Not archived"
            )
        except ClientError as e:
            print(f"Error checking thaw status for {key}: {e}")

    def thaw(self, bucket, key, tier="Standard", days=7):
        try:
            self.client.restore_object(
                Bucket=bucket,
                Key=key,
                RestoreRequest={"Days": days, "GlacierJobParameters": {"Tier": tier}},
            )
            print(f"Restore request initiated for {key} (Tier: {tier}, Days: {days})")
        except ClientError as e:
            print(f"Error restoring {key}: {e}")
