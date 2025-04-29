# PYTHON_ARGCOMPLETE_OK
import sys
import argparse
import argcomplete

# import project command-line interfaces
from a_series_of_tubes.ec2manager import cli as ec2cli


def main():
    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatics Meta-Framework (cbmf)"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    ec2manager_parser = subparsers.add_parser("ec2", help="Manage EC2 instances")
    ec2manager_parser.set_defaults(func=dispatch_ec2manager)

    # enable tab completion
    argcomplete.autocomplete(parser)

    args, unknown = parser.parse_known_args()
    args.func(unknown)


def dispatch_ec2manager(argv):
    sys.argv = [sys.argv[0]] + argv
    ec2cli.main()
