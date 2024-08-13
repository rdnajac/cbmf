import sys
import subprocess
from typing import List, Optional, Callable, Dict
from pathlib import Path

from .utils.cli import create_parser
from .utils.colorprinter import ColorPrinter as pr
from .utils.genomemanager import download_genome_file, download_all_files_for_species
# call this with human or mouse based on args
from .tests.test_colorprinter import smoke_test

def run_quality_control(args):
    pr.info(f"Running quality control on {args.input_directory} and saving results to {args.output_directory}.")
    # TODO: Implement actual QC logic

def run_alignment(args):
    pr.info(f"Running {args.aligner} alignment for {args.species} genome.")
    pr.info(f"Input directory: {args.input_directory}")
    pr.info(f"Output directory: {args.output_directory}")
    # TODO: Implement actual alignment logic

def initialize_pipeline(args):
    pr.info("Initializing pipeline files for {args.species}")
    download_all_files_for_species(args.species)

def check_pipeline_status(args):
    pr.info("Checking pipeline status...")
    # TODO: Implement status check

def run_test_suite(args):
    pr.info("Running test suite")
    smoke_test()

def main(argv: Optional[List[str]] = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)

    command_handlers = {
        'download': lambda args: download_genome_file(args.species, args.file),
        'init': initialize_pipeline,
        'qc': run_quality_control,
        'align': run_alignment,
        'status': check_pipeline_status,
        'test': run_test_suite
    }

    try:
        handler = command_handlers.get(args.command)
        if handler:
            pr.info(f"Executing {args.command}:")
            handler(args)
            pr.success(f"{args.command} completed successfully.")
        else:
            pr.error(f"Unknown command: {args.command}")
            return 1
    except ValueError as e:
        pr.error(f"Error: {str(e)}")
        return 1
    except subprocess.CalledProcessError as e:
        pr.error(f"Command failed: {str(e)}")
        return 1
    except Exception as e:
        pr.error(f"An unexpected error occurred: {str(e)}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
