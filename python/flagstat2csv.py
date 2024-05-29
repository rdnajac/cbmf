#! /usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from colorprinter import ColorPrinter as pr


def bail(msg):
    """Print an error message and exit the program"""
    pr.error(msg)
    sys.exit(1)


# Define regular expressions to match the flagstat output
PATTERNS = {
    "total_reads": re.compile(r"(\d+) \+ (\d+) in total"),
    "mapped": re.compile(r"(\d+) \+ (\d+) mapped"),
    "duplicates": re.compile(r"(\d+) \+ (\d+) duplicates"),
}


def parse_stats(data):
    """
    Parse the flagstat output and return a dictionary with the results

    Args:
        data (str): The flagstat output

    Returns:
        dict: A dictionary with the parsed results
    """
    results = {}
    try:
        for key, pattern in PATTERNS.items():
            match = pattern.search(data)
            if match:
                results[key] = {
                    "qc_passed": int(match.group(1)),
                    "qc_failed": int(match.group(2)),
                }
        return results
    except Exception as e:
        bail(f"Error parsing flagstat output: {e}")


def to_percent(part, whole):
    """Calculate the percentage of part in whole"""
    return (part / whole) * 100 if whole else 0


def calculate_additional_stats(stats, total_reads):
    """
    Calculate additional statistics based on the parsed flagstat output.

    Args:
        stats (dict): The parsed flagstat output
        total_reads (int): The total number of reads

    Returns:
        dict: An updated dictionary with the additional statistics
    """
    mapped_reads = stats.get("mapped", {}).get("qc_passed", 0)
    duplicate_reads = stats.get("duplicates", {}).get("qc_passed", 0)
    unique_reads = mapped_reads - duplicate_reads
    return {
        "% Mapped Reads": to_percent(mapped_reads, total_reads),
        "% Duplicate Reads": to_percent(duplicate_reads, total_reads),
        "Unique Reads": unique_reads,
        "% Unique Reads": to_percent(unique_reads, total_reads),
        "% Unique of Mapped": to_percent(unique_reads, mapped_reads),
    }


def process_bam_file(file):
    """
    Process a BAM file and return the statistics by running
    samtools subprocesses to get the flagstats and total read count.

    Args:
        file (str): The path to the BAM file

    Returns:
        tuple: A tuple with the filename and the statistics
    """
    try:
        pr.info(f"Processing {file}")
        subprocess.run("samtools --version", shell=True, check=True)

        flagstat_cmd = subprocess.run(
            f"samtools flagstat -@ 32 {file}",
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )
        total_reads_cmd = subprocess.run(
            f"samtools view -@ 32 -c {file}",
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )

        stats = parse_stats(flagstat_cmd.stdout)
        total_reads = int(total_reads_cmd.stdout.strip())
        additional_stats = calculate_additional_stats(stats, total_reads)
        stats.update(additional_stats)

        return os.path.basename(file), stats

    except subprocess.CalledProcessError as e:
        pr.error(f"Error running samtools: {e}")
        raise RuntimeError(
            f"Failed to process file {file} due to samtools error."
        ) from e

    return os.path.basename(file), stats


def collect_bam_stats(directory, nproc=4):
    bam_files = glob.glob(os.path.join(directory, "*.bam"))
    results_dict = {}

    with ThreadPoolExecutor(max_workers=nproc) as executor:
        futures = {executor.submit(process_bam_file, file): file for file in bam_files}
        for future in as_completed(futures):
            filename, stats = future.result()
            results_dict[filename] = stats

    return results_dict


HEADERS = [
    "Sample ID",
    "Filename",
    "Total Reads",
    "Passing Reads",
    "% Passed Reads",
    "Mapped Reads",
    "% Mapped Reads",
    "Duplicate Reads",
    "% Duplicate Reads",
    "Unique Reads",
    "% Unique Reads",
    "% Unique of Mapped",
]


def write_csv(results_dict, filename="flagstat.csv"):
    """
    Write the results to a CSV file.

    Args:
        results_dict (dict): A dictionary containing the results
        filename (str): The name of the output CSV file
    """

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(HEADERS)
        for filename, stats in results_dict.items():
            total_reads = stats.get("total_reads", {}).get("qc_passed", 0) + stats.get(
                "total_reads", {}
            ).get("qc_failed", 0)
            row = [
                filename.split(".")[0],
                filename,
                total_reads,
                stats.get("total_reads", {}).get("qc_passed", 0),
                to_percent(
                    stats.get("total_reads", {}).get("qc_passed", 0), total_reads
                ),
                stats.get("mapped", {}).get("qc_passed", 0),
                stats.get("% Mapped Reads", 0),
                stats.get("duplicates", {}).get("qc_passed", 0),
                stats.get("% Duplicate Reads", 0),
                stats.get("Unique Reads", 0),
                stats.get("% Unique Reads", 0),
                stats.get("% Unique of Mapped", 0),
            ]
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(
        description="Parse flagstat output and write the results to a CSV file."
    )
    parser.add_argument(
        "directory",
        type=str,
        help="The directory containing the BAM files",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="flagstat.csv",
        help="The name of the output CSV file",
    )
    parser.add_argument(
        "-n",
        "--nproc",
        type=int,
        default=32,
        help="The number of processes to use for parallel processing",
    )
    args = parser.parse_args()

    try:
        results_dict = collect_bam_stats(args.directory, args.nproc)
        write_csv(results_dict, args.output)
        pr.success(f"Results written to {args.output}")
    except Exception as e:
        pr.error(f"Error processing BAM files: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
