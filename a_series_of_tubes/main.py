import sys
import subprocess
import importlib
from typing import List, Optional, Callable, Dict
from pathlib import Path

from .utils.cli import create_parser
from .utils.colorprinter import ColorPrinter as pr
from .utils.genomemanager import (
    download_genome_file,
    download_all_files_for_species,
    decompress_files,
)
from .utils.run_script import run_script
from .utils.progressbar import ProgressBar

from .tests.test_colorprinter import smoke_test
from .tests.test_progressbar import *

from .aligners import align_reads


def run_quality_control(args):
    pr.info(f" QC in: {args.input_directory}, QC out: {args.output_directory}")
    # TODO: Implement actual QC logic


# dict mapping human/mouse to full release name
REFERENCE = {"human": "GRCh38", "mouse": "GRCm38"}


def run_alignment(args):
    pr.info(f"Running {args.aligner} alignment for {args.species} genome.")
    pr.info(f"Input directory: {args.input_directory}")
    pr.info(f"Output directory: {args.output_directory}")

    # Get the list of FASTQ files in the input directory
    # fastq_files = sorted(Path(args.input_directory).glob("*.fastq.gz"))
    # if len(fastq_files) < 2:
    #     pr.error("Not enough FASTQ files found. Need at least a pair of files.")
    #     return

    # Assume the first two files are the paired-end reads
    r1_fastq, r2_fastq = ["r1.fastq.gz", "r2.fastq.gz"]
    reference = "ref"
    should_be_none = align_reads(args.aligner, r1_fastq, r2_fastq, reference)
    if should_be_none is not None:
        print('huh?')

        # Save the aligned output to a file
        # output_sam = Path(args.output_directory) / f"{args.species}_{args.aligner}_aligned.sam"
        # with output_sam.open('w') as f:
        #     f.write(aligned_output)

        # pr.success(f"Alignment completed. Output saved to {output_sam}")


def initialize_pipeline(args):
    spec_genome_dir = Path.home() / "cbmf" / "genomes" / args.species
    spec_genome_dir.mkdir(parents=True, exist_ok=True)
    pr.info(f"Downloading genome files for {args.species} to {spec_genome_dir}")
    download_all_files_for_species(args.species)
    decompress_files(spec_genome_dir)


def check_pipeline_status(args):
    pr.info("Checking pipeline status...")
    # TODO: Implement status check


def run_test_suite(args):
    if not args.tests:
        smoke_test()
        # print(args)
        # pr.info("Running all tests")
        # Implement logic to run all tests
    else:
        # test_progressbar()
        # test_two_progressbars()
        exit()
        pr.info(f"Running specified tests: {', '.join(args.tests)}")
        for test in args.tests:
            run_single_test(test)


def run_single_test(test_name):
    test_file = f"test_{test_name}.py"
    test_path = Path(__file__).parent / "tests" / test_file

    if not test_path.exists():
        pr.error(f"Test file {test_file} not found in path {test_path}")
        return

    try:
        module_name = f"tests.{test_file[:-3]}"
        module = importlib.import_module(module_name)
        if hasattr(module, "main"):
            module.main()
        else:
            pr.error(f"No main function found in {test_file}")
    except ImportError as e:
        pr.error(f"Error importing test module: {str(e)}")
    except Exception as e:
        pr.error(f"Error running test: {str(e)}")


def main(argv: Optional[List[str]] = None) -> int:

    # if no args are passed, print help
    if len(sys.argv) == 1:
        run_script("~/cbmf/scripts/SUCCESS")
        return 0

    parser = create_parser()
    args = parser.parse_args(argv)

    command_handlers = {
        "download": lambda args: download_genome_file(args.species, args.file),
        "init": initialize_pipeline,
        "qc": run_quality_control,
        "align": run_alignment,
        "status": check_pipeline_status,
        "test": run_test_suite,
    }

    try:
        handler = command_handlers.get(args.command)
        if handler:
            pr.warning(f"{args}")
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
