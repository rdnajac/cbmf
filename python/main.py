#!/usr/bin/env python3

import argparse
from typing import Dict, Set


# Define mappings for valid species and aligners
VALID_SPECIES_MAPPING: Dict[str, Set[str]] = {
    "human": {"hu", "h"},
    "mouse": {"mo", "m"},
}

VALID_ALIGNERS_MAPPING: Dict[str, Set[str]] = {
    "bwa": {"bwa"},
    "hisat2": {"hisat2", "hisat"},
    "bowtie2": {"bowtie2", "bowtie"},
    "star": {"star"},
    "subread-align": {"subread-align", "subread"},
    "subjunct": {"subjunct"},
}


def validate(value: str, mappings: Dict[str, Set[str]]) -> str:
    """
    Validate and normalize a string value based on a mapping of valid options.

    This function converts the input value to lowercase and checks it against
    the provided mappings.  If the value matches one of the aliases in the
    mappings, it returns the canonical name.  If not, it raises a ValueError.

    Args:
        value (str): The value to validate.
        mappings (Dict[str, Set[str]]): A dictionary where keys are canonical
            names and values are sets of aliases.

    Returns:
        str: The canonical name corresponding to the input value.

    Raises:
        ValueError: If the value does not match any of the valid options in the mappings.
    """
    value = value.lower()
    for canonical_name, aliases in mappings.items():
        if value in aliases:
            return canonical_name
    raise ValueError("Value not supported")


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments and validate them.

    This function uses argparse to parse command-line arguments, validate the 
    species and aligner arguments, and ensure there are at least two additional 
    arguments provided. If fewer than two additional arguments are given, it 
    extends the list with `None` values to ensure the correct number of arguments.

    Returns:
        argparse.Namespace: An object containing the parsed arguments and their values.
    """
    parser = argparse.ArgumentParser(description="Python application template")
    parser.add_argument(
        "-s", "--species", type=str, help="Species (e.g., human, mouse)"
    )
    parser.add_argument(
        "-a",
        "--aligner",
        type=str,
        help="Aligner (e.g., bwa, bowtie2, star, subread-align, subjunct)",
    )
    parser.add_argument("args", nargs="*", help="Additional arguments")

    args = parser.parse_args()

    # Validate and normalize species
    if args.species:
        args.species = validate(args.species, VALID_SPECIES_MAPPING)

    # Validate and normalize aligner
    if args.aligner:
        args.aligner = validate(args.aligner, VALID_ALIGNERS_MAPPING)

    # Ensure there are at least two additional arguments
    if len(args.args) < 2:
        args.args.extend([None] * (2 - len(args.args)))

    return args


def print_args(args: argparse.Namespace) -> None:
    """
    Print the parsed arguments in a formatted manner.

    This function prints the values of the species, aligner, and additional arguments.
    It ensures that all provided arguments are displayed clearly.

    Args:
        args (argparse.Namespace): An object containing the parsed arguments and their values.
    """
    print(f"Species: {args.species if args.species else 'Not specified'}")
    print(f"Aligner: {args.aligner if args.aligner else 'Not specified'}")
    print(f"Additional Argument 1: {args.args[0]}")
    print(f"Additional Argument 2: {args.args[1]}")


def main():
    """Main function to handle command-line arguments and execute the script."""
    args = parse_args()
    print_args(args)


if __name__ == "__main__":
    main()
