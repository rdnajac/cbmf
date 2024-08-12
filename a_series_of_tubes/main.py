import sys
import argparse
from typing import List, Optional, Callable, Dict
from pathlib import Path
from .utils.colorprinter import ColorPrinter as pr
from .utils.genomemanager import download_genome_file, PIPELINE_FILES
from .tests.test_colorprinter import smoke_test

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='cbmf', description="Combinatorial Bioinformatic Meta-Framework (CBMF)")
    
    # Species selection
    parser.add_argument(
        '--species', 
        choices=['human', 'mouse'], 
        help='Species to use (human or mouse)'
    )
    parser.add_argument(
        '--hu', '--human', 
        action='store_const', 
        dest='species', 
        const='human', 
        help='Shortcut for human species'
    )
    parser.add_argument(
        '--mo', '--mouse', 
        action='store_const', 
        dest='species', 
        const='mouse', 
        help='Shortcut for mouse species'
    )

    # Parent parser for common arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-i", "--input-directory", type=Path, help="Input directory path")
    parent_parser.add_argument("-o", "--output-directory", type=Path, help="Output directory path")

    # Subparsers
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Align subcommand
    parser_align = subparsers.add_parser('align', parents=[parent_parser], help='Run alignment')
    align_method = parser_align.add_mutually_exclusive_group(required=True)
    align_method.add_argument('--bwa', action='store_true', help='Use BWA aligner')
    align_method.add_argument('--hisat2', action='store_true', help='Use HISAT2 aligner')
    align_method.add_argument('--bowtie', action='store_true', help='Use Bowtie aligner')

    # QC subcommand
    subparsers.add_parser('qc', parents=[parent_parser], help='Run quality control')

    # Test subcommand
    subparsers.add_parser('test', help='Run test suite')

    # Download subcommand
    parser_download = subparsers.add_parser('download', help='Download genome files')
    parser_download.add_argument('--file', required=True, help='File to download')

    return parser

def download_helper(args: argparse.Namespace) -> None:
    try:
        species = args.species
        if not species:
            raise ValueError("Species not specified. Use --species, --hu, or --human for human, or --mo, or --mouse for mouse.")
        pr.info(f"Initializing download for {args.file} for {species}.")
        download_genome_file(species, args.file)
        pr.success(f"{args.file} download completed.")
    except ValueError as e:
        pr.error(f"Error: {str(e)}")
        raise
    except Exception as e:
        pr.error(f"An unexpected error occurred: {str(e)}")
        raise

def run_quality_control(args: argparse.Namespace) -> None:
    pr.info(f"Running quality control on {args.input_directory} and saving results to {args.output_directory}.")
    # Implement the actual QC logic here

def run_alignment(args: argparse.Namespace) -> None:
    pr.info(f"Aligning {args.input_directory} to {args.reference} and saving results to {args.output_directory}.")
    
    # Dictionary of aligners
    aligners: Dict[str, Callable[[Path, Path, Path], None]] = {
        'bwa': "align_with_bwa",
        'hisat2': "align_with_hisat2",
        'bowtie': "align_with_bowtie"
    }
    pr.info(f"Using {args.command} aligner.")

def run_test_suite(args: argparse.Namespace) -> None:
    pr.info("Running test suite")
    smoke_test()

COMMAND_HANDLERS: Dict[str, Callable[[argparse.Namespace], None]] = {
    "download": download_helper,
    "qc": run_quality_control,
    "align": run_alignment,
    "test": run_test_suite,
}

def main(argv: Optional[List[str]] = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    
    if not args.command:
        parser.print_help()
        return 0

    try:
        handler = COMMAND_HANDLERS.get(args.command)
        if handler:
            handler(args)
        else:
            raise ValueError(f"Unknown command: {args.command}")
    except Exception as e:
        pr.error(str(e))
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
