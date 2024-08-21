import argparse
from pathlib import Path


def create_parser() -> argparse.ArgumentParser:
    """Create the main CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatic Meta-Framework (CBMF)"
    )

    # The species_parser lets the user specify the species with aliases
    species_parser = argparse.ArgumentParser(add_help=False)
    # Species arg is mutually exclusive so it fails if more than one is provided
    # species_group = species_parser.add_mutually_exclusive_group(required=True)
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
    # Add more species as needed...

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
        choices=[
            "bowtie_index",
            "bwa_index",
            "samtools_index",
            "fasta",
            "hisat2_index",
            "refseq_gff",
            "refseq_gtf",
        ],
        help="choose file(s) to download",
    )

    # Align reads to a reference genome
    subparser_align = subparsers.add_parser(
        "align",
        help="run alignment",
        # "align", parents=[species_parser], help="run alignment"
    )
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
        # TODO attempt to resolve from the species and aligner
        help="Path to the reference genome index",
    )

    test_parser = subparsers.add_parser("test", help="run test suite")
    test_parser.add_argument("tests", nargs="*", help="list of tests to run")

    return parser


def validate_args(args):
    """
    Validate parsed arguments.

    First, check if an arg exists, then check it against a set of conditions.
    Raise erros here and not everywhere else
    Try to centralize error handling.
    """

    # The input directory must exist and be a directory
    if args.input_dir.exists():
        if not args.input_dir.is_dir():
            raise NotADirectoryError(f"Input directory not found: {args.input_dir}")
        # It can't be empty either
        if not list(args.input_dir.iterdir()):
            raise ValueError(f"Input directory is empty: {args.input_dir}")

    # If the output directory doesn't exist, create it
    if not args.output_dir.exists():
        args.output_dir.mkdir(parents=True, exist_ok=True)


def parse_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = create_parser()
    args = parser.parse_args(argv)
    validate_args(args)
    return args
