import sys
import argparse
from .utils.colorprinter import ColorPrinter as pr
from .utils.genomemanager import download_genome_file
from .tests.test_colorprinter import smoke_test

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Combinatorial Bioinformatic Meta-Framework (CBMF) Command-Line Interface",
    )

    parser.add_argument("-V", "--version", action="version", version="CBMF v0.9")


    command = parser.add_mutually_exclusive_group(required=True)
    command.add_argument("--init", action="store_true")
    command.add_argument("--qc", action="store_true")
    command.add_argument("--test", action="store_true")
    command.add_argument("-d", "--download", choices=["bwa_index", "samtools_index", "fasta", "hisat2_index", "refseq_gff", "refseq_gtf"])

    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-r", "--reference")

    species = parser.add_mutually_exclusive_group()
    species.add_argument("--human", action="store_true")
    species.add_argument("--mouse", action="store_true")

    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    group.add_argument("-q", "--quiet", action="store_true", help="suppress all output")

    return parser.parse_args()

def download_helper(options):
    try:
        if options.human:
            species = "human"
        elif options.mouse:
            species = "mouse"
        else:
            raise ValueError("Species not specified. Use --human or --mouse.")
        pr.info(f"Initializing {options.download} for {species}.")
        download_genome_file(species, options.download)
        pr.success(f"{options.download} initialization completed.")
    except ValueError as e:
        raise ValueError(f"Error: {str(e)}")
    except Exception as e:
        raise Exception(f"An unexpected error occurred: {str(e)}")

def run_quality_control(options):
    pr.info(
        f"Running quality control on {options.input} and saving results to {options.output}."
    )

def run_test_suite():
    smoke_test()

def main():
    try:
        args = parse_arguments()
        
        if args.download:
            download_helper(args)
        elif args.qc:
            run_quality_control(args)
        elif args.test:
            run_test_suite()

    except ValueError as e:
        pr.error(f"Error: {str(e)}")
        return 1
    except Exception as e:
        pr.error(f"An unexpected error occurred: {str(e)}")
        return 2

    return 0

if __name__ == "__main__":
    sys.exit(main())
