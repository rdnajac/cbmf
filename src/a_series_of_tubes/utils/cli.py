import argparse
import logging
from pathlib import Path
from ..config import FILES


def create_parser() -> argparse.ArgumentParser:
    """Create the main CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatic Meta-Framework (CBMF)"
    )

    # Add verbosity option
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase output verbosity (e.g., -v, -vv, -vvv)",
    )

    # The species_parser lets the user specify the species with aliases
    species_parser = argparse.ArgumentParser(add_help=False)
    species_group = species_parser.add_mutually_exclusive_group(required=True)
    species_group.add_argument(
        "--species",
        choices=["human", "mouse"],
        help="choose a species",
    )
    species_group.add_argument(
        "--hu",
        "--human",
        action="store_const",
        dest="species",
        const="human",
        help='shorthand for "--species human"',
    )
    species_group.add_argument(
        "--mo",
        "--mouse",
        action="store_const",
        dest="species",
        const="mouse",
        help='shorthand for "--species mouse"',
    )

    # Create subparsers for each command and add them to the main parser
    subparsers = parser.add_subparsers(dest="command", required=True)

# Download a specific genome file from the NCBI seqs for pipelines
    subparser_download = subparsers.add_parser(
        "download",
        aliases=["dl"],
        parents=[species_parser],
        help="download NCBI genome files",
    )
    subparser_download.add_argument(
        "files",
        nargs="+",
        choices=list(FILES.keys()) + ["ALL"],
        help="choose file(s) to download",
    )

    # Align reads to a reference genome
    subparser_align = subparsers.add_parser("align", help="run alignment")
    subparser_align.add_argument(
        "aligner", choices=["bwa", "hisat2", "bowtie2"], help="choose aligner"
    )
    subparser_align.add_argument(
        "-i",
        "--input-directory",
        type=lambda p: Path(p).expanduser().resolve(),
        required=True,
        help="Input directory containing FASTQ files",
    )
    subparser_align.add_argument(
        "-o",
        "--output-directory",
        type=lambda p: Path(p).expanduser().resolve(),
        default=Path.cwd(),
        help="Output directory for SAM files (default: current directory)",
    )
    subparser_align.add_argument(
        "-r",
        "--reference",
        type=lambda p: Path(p).expanduser().resolve(),
        help="Path to the reference genome index",
    )

    # Test command
    subparsers.add_parser("test", help="run test suite")

    return parser


def validate_args(args):
    """
    Validate parsed arguments.
    """
    if hasattr(args, "input_directory") and args.input_directory:
        if not args.input_directory.exists():
            raise NotADirectoryError(
                f"Input directory not found: {args.input_directory}"
            )
        if not args.input_directory.is_dir():
            raise NotADirectoryError(
                f"Input path is not a directory: {args.input_directory}"
            )
        if not list(args.input_directory.iterdir()):
            raise ValueError(f"Input directory is empty: {args.input_directory}")

    if hasattr(args, "output_directory") and args.output_directory:
        if not args.output_directory.exists():
            args.output_directory.mkdir(parents=True, exist_ok=True)


def parse_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = create_parser()
    args = parser.parse_args(argv)
    validate_args(args)

    # Set log level based on verbosity
    verbosity_map = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    args.log_level = verbosity_map.get(args.verbose, logging.DEBUG)

    return args
