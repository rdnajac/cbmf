import sys
import subprocess
import importlib
from pathlib import Path

from .aligners import align_reads

from .utils import parse_args, pr, run_script, download_helper


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]  # Skip the script name

    if len(argv) == 0:
        run_script("~/cbmf/scripts/SUCCESS")
        return 0

    args = parse_args(argv)
    # print(args)

    subcommands = {
        "download": lambda args: download_helper(args.species, args.files),
        # "init": initialize_pipeline,
        # "qc": run_quality_control,
        # "align": run_alignment,
        # "status": check_pipeline_status,
        # "test": run_test_suite,
    }

    try:
        handler = subcommands.get(args.command)
        if handler:
            pr.warning(f"{args}")
            handler(args)
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
        print("huh?")

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


if __name__ == "__main__":
    sys.exit(main())
