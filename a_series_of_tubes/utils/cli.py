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

    # directory_parser = argparse.ArgumentParser(add_help=False)
    # directory_parser.add_argument(
    #     "-i", "--in",
    #     "--input-directory", type=Path, metavar="DIR", required=True
    # )
    # directory_parser.add_argument(
    #     "-o", "--out" "--output-directory", type=Path, metavar="DIR", default=Path(".")
    # )

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

    subparser_align = subparsers.add_parser(
        "align", parents=[species_parser], help="run alignment"
    )
    subparser_align.add_argument(
        "aligner", choices=["bwa", "hisat2", "bowtie"], help="choose aligner"
    )

    # subparsers.add_parser(
    #     "init", parents=[species_parser], help="initialize pipeline ref/idx files"
    # )

    # subparsers.add_parser("qc", parents=[directory_parser], help="run quality control")
    # subparsers.add_parser("status", help="check pipeline status")

    # test_parser = subparsers.add_parser("test", help="run test suite")
    # test_parser.add_argument("tests", nargs="*", help="list of tests to run")

    return parser


def parse_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = create_parser()
    return parser.parse_args(argv)
