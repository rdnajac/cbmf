# utils/cli.py
import argparse
from typing import Dict, Any


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
    # parser.add_argument("--test",  action="store_true", help="Run the test suite")

    # # Common arguments
    # parser.add_argument("-i", "--input",     help="Input file/directory")
    # parser.add_argument("-o", "--output",    help="Output directory")
    # parser.add_argument("-r", "--reference", help="Path to reference file")
    # parser.add_argument("-s", "--species", choices=["human", "mouse"])

    # # Alignment-specific arguments
    # parser.add_argument("--aligner", choices=["bwa", "star", "hisat2", "bowtie2"])

    return parser.parse_args()


def validate_options(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Validate and process the parsed command-line arguments.

    This function checks for required arguments, mutually exclusive options,
    and performs any necessary preprocessing of the arguments.

    :param args: The parsed command-line arguments
    :type args: argparse.Namespace
    :return: A dictionary of validated and processed options
    :rtype: Dict[str, Any]
    :raises ValueError: If the provided arguments are invalid or inconsistent
    """
    pass
    options = vars(args)
    commands = ["init", "qc", "align", "test"]
    active_commands = [cmd for cmd in commands if options[cmd]]

    if len(active_commands) == 0:
        raise ValueError(
            "Please specify a command or rerun with --help for more information."
        )

    elif len(active_commands) > 1:
        raise ValueError("Multiple commands specified. Please choose only one.")

    #     if options["init"]:
    #         if not options["species"]:
    #             raise ValueError("--init command requires --species argument.")

    #     elif options["qc"]:
    #         if not options["input"] or not options["output"]:
    #             raise ValueError(
    #                 "--qc command requires both --input and --output arguments."
    #             )

    #         elif options["align"]:
    #         if not options["input"] or not options["output"] or not options["reference"] or not options["aligner"]:
    #             raise ValueError(
    #                 "--align command requires --input, --output, --reference, and --aligner arguments."
    #             )

    return options
