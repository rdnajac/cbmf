# tests/test_cli.py
from ..utils.cli import parse_arguments 

def print_parsed_args():
    args = parse_arguments()
    print(args)
