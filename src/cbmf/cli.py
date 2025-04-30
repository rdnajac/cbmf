#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argparse
import argcomplete

from a_series_of_tubes.ec2manager import cli as ec2cli
from a_series_of_tubes.s3manager import cli as s3cli


dispatch = {
    "ec2": ec2cli.main,
    "s3": s3cli.main,
}


def main():
    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatics Meta-Framework (cbmf)"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, func in dispatch.items():
        sp = subparsers.add_parser(name, add_help=False)  # suppress dummy -h
        sp.set_defaults(func=func)

    argcomplete.autocomplete(parser)
    args, unknown = parser.parse_known_args()
    if "-h" in unknown or "--help" in unknown:
        unknown.append("--help")  # ensure forwarded
    args.func(unknown)


if __name__ == "__main__":
    main()
