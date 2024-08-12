import sys
from .utils.cli import parse_arguments, validate_options
from .utils.colorprinter import ColorPrinter as pr
# from .utils.genomemanager import download_genome_file
from .tests.test_colorprinter import smoke_test


def run_quality_control(options):
    pr.info(
        f"Running quality control on {options['input']} and saving results to {options['output']}."
    )
    # TODO: Implement actual QC functionality here


def run_alignment(options):
    pr.info(
        f"Aligning {options['input']} to {options['reference']} using {options['aligner']} and saving results to {options['output']}."
    )
    # TODO: Implement actual alignment functionality here


def download_genome_files(options):
    pr.info(f"Initializing genome files for {options['species']}.")
    # TODO: Implement genome file download logic here


def main():
    try:
        args = parse_arguments()
        # options = validate_options(args)
        print(args)

#         if options["init"]:
#             download_genome_files(options)

#         elif options["qc"]:
#             run_quality_control(options)

#         elif options["align"]:
#             run_alignment(options)

#         elif options["test"]:
#             smoke_test()

    except ValueError as e:
        pr.error(f"Error: {str(e)}")
        return 1
    except Exception as e:
        pr.error(f"An unexpected error occurred: {str(e)}")
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
