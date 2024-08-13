import argparse
from pathlib import Path


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatic Meta-Framework (CBMF)"
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-i", "--input-directory", type=Path, metavar="DIR")
    parent_parser.add_argument("-o", "--output-directory", type=Path, metavar="DIR")

    species_parser = argparse.ArgumentParser(add_help=False)
    species_group = species_parser.add_mutually_exclusive_group()
    species_group.add_argument("--species", choices=["human", "mouse"])
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

    subparsers = parser.add_subparsers(dest="command", required=True)

    subparser_download = subparsers.add_parser(
        "download", parents=[species_parser], help="download genome files"
    )
    subparser_download.add_argument(
        "file",
        choices=[
            "bowtie_index",
            "bwa_index",
            "samtools_index",
            "fasta",
            "hisat2_index",
            "refseq_gff",
            "refseq_gtf",
        ],
        help="choose file to download",
    )

    subparsers.add_parser(
        "init", parents=[species_parser], help="initialize pipeline ref/idx files"
    )

    subparser_align = subparsers.add_parser(
        "align", parents=[parent_parser, species_parser], help="run alignment"
    )
    subparser_align.add_argument(
        "aligner", choices=["bwa", "hisat2", "bowtie", "star"], help="choose aligner"
    )

    subparsers.add_parser("qc", parents=[parent_parser], help="run quality control")
    subparsers.add_parser("status", help="check pipeline status")

    subparsers.add_parser("test", help="run test suite")
    return parser
