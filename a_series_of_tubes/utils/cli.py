# utils/cli.py
import argparse
# from typing import Dict, Any


def parse_arguments():
    """
    Parse command-line arguments for the frontend script.

    This function sets up the argument parser and defines all the command-line
    options for the Combinatorial Bioinformatic Meta-Framework (CBMF) tool.

    :return: A namespace object containing the parsed arguments
    :rtype: argparse.Namespace
    :raises SystemExit: If the user provides invalid arguments
    """

    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatic Meta-Framework (CBMF) Command-Line Interface",
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # General arguments
    parser.add_argument("-V", "--version", action="version", version="CBMF v0.9")

    # Verbosity options
    # define a group of mutually exclusive options
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")

    # # Command selection
    # parser.add_argument("--init",  action="store_true", help="Initialize files")
    # parser.add_argument("--qc",    action="store_true", help="Run quality control")
    # parser.add_argument("--align", action="store_true", help="Align reads to reference")
    parser.add_argument("--test",  action="store_true", help="Run the test suite")

    # # Common arguments
    parser.add_argument("-i", "--input",     help="Input file/directory")
    parser.add_argument("-o", "--output",    help="Output directory")
    parser.add_argument("-r", "--reference", help="Path to reference")
    
    species = parser.add_mutually_exclusive_group()
    species.add_argument("--human", action="store_true")
    species.add_argument("--mouse", action="store_true")

    # # Alignment-specific arguments
    # parser.add_argument("--aligner", choices=["bwa", "star", "hisat2", "bowtie2"])

    return parser.parse_args()

