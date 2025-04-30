#!/usr/bin/env python3

import argparse
from a_series_of_tubes.s3manager.s3_helper import S3Helper


dispatch = {
    "list": lambda s3, b, p, **_: s3.list_files(b, p),
    "metadata": lambda s3, b, p, **_: [
        print(s3.get_metadata(b, k)) for k in s3.list_files(b, p)
    ],
    "thaw-status": lambda s3, b, p, **_: [
        print(f"{k} - Storage Class: {status[0]}, Restore Status: {status[1]}")
        for k in s3.list_files(b, p)
        if (status := s3.get_thaw_status(b, k))
    ],
    "thaw": lambda s3, b, p, **kw: [
        s3.thaw(b, k, tier=kw["tier"], days=kw["days"]) for k in s3.list_files(b, p)
    ],
}


def main(argv=None):
    parser = argparse.ArgumentParser(description="S3 Glacier Management Tool")
    parser.add_argument("uri", help="S3 URI (e.g., s3://bucket/prefix/)")
    parser.add_argument("action", choices=dispatch.keys(), help="Action to perform")
    parser.add_argument(
        "--tier", default="Standard", help="Thaw tier (Expedited/Standard/Bulk)"
    )
    parser.add_argument(
        "--days", type=int, default=7, help="Days to keep thawed files available"
    )

    args = parser.parse_args(argv)
    s3 = S3Helper()
    bucket, prefix = s3.parse_uri(args.uri)
    dispatch[args.action](s3, bucket, prefix, tier=args.tier, days=args.days)


if __name__ == "__main__":
    main()
