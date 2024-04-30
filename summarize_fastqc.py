#!/usr/bin/env python3
""" 
Summarize all the FastQC reports in a directory to a CSV file.
Works with FastQC from .fastq.gz files and .bam files.

Usage: 
    python summarize_fastqc.py
"""

import os
import csv
import re
import zipfile

summary_fields = [
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
]


def summarize_zipped_reports(directory):
    """
    Summarize all the FastQC reports in a directory to a CSV file.

    Args:
    directory: the directory containing the FastQC zip files

    Returns:
    A list of dictionaries containing the summary fields and status.
    """
    reports = []
    for root, dirs, files in os.walk(directory):
        for fname in files:
            if fname.endswith("fastqc.zip"):
                try:
                    with zipfile.ZipFile(os.path.join(root, fname), "r") as zip_ref:
                        summary_content = next(
                            (
                                zip_ref.read(file).decode("utf-8")
                                for file in zip_ref.namelist()
                                if "summary.txt" in file
                            ),
                            None,
                        )
                        if summary_content:
                            qc_summary = summarize(summary_content.splitlines(), fname)
                            reports.append(qc_summary)
                except zipfile.BadZipFile:
                    print(f"Bad zip file: {fname}")
                except Exception as e:
                    print(f"Failed to process {fname}: {e}")
    return reports


def summarize(summary_lines, fname):
    """
    Summarize the FastQC report summary content.

    Args:
    summary_lines: list of lines from the FastQC summary.txt file
    fname: the filename of the FastQC report

    Returns:
    A dictionary with the summary fields as keys and the status as values.
    """
    summary_dict = {line.split("\t")[1]: line.split("\t")[0] for line in summary_lines}
    summary_dict["Sample"] = re.sub(r"_001\.fastqc\.zip", "", fname)
    return summary_dict


def write_summary_to_csv(reports, output_file):
    try:
        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Sample"] + summary_fields)
            for report in reports:
                status_list = [report.get(field, "") for field in summary_fields]
                writer.writerow([report["Sample"]] + status_list)
    except Exception as e:
        print(f"Error writing to CSV: {e}")


import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Summarize all the FastQC reports in a directory to a CSV file."
    )
    parser.add_argument(
        "directory",
        type=str,
        help="The directory containing the FastQC zip files",
        default=os.getcwd(),
        nargs="?",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="The output CSV file",
        default="fastqc_summary.csv",
        nargs="?",
    )
    args = parser.parse_args()
    reports = summarize_zipped_reports(args.directory)
    write_summary_to_csv(reports, args.output_file)


"""
to run with default output file name
python summarize_fastqc.py <directory>
./summarize_fastqc.py <directory>
"""
if __name__ == "__main__":
    main()
