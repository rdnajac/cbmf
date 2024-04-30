#!/usr/bin/env python3

import os
import csv
import sys
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

data_fields = [
    "file",
    "File type",
    "Encoding",
    "Total Sequences",
    "Total Bases",
    "% Bases >= 30",
    "Sequences flagged as poor quality",
    "Mean Quality Score",
    "Sequence length",
    "Mapped Reads",
    "% Mapped Reads",
    "Duplicate Reads",
    "% Duplicate Reads (of all reads)",
    "Unique Reads",
    "% Unique Reads (of all reads)",
    "% Unique (of mapped)",
]


class QCManager:
    """Class to manage and summarize FastQC reports."""

    def __init__(self, directory):
        self.summaries = {}
        self.data_dict = {}
        self.process_zipped_reports(directory)

    def summarize(self, summary_txt):
        """Parse `summary.txt` file and store the results in a dictionary."""
        for line in summary_txt:
            status, field, filename = line.decode("utf-8").strip().split("\t")
            if filename not in self.summaries:
                self.summaries[filename] = {}
            self.summaries[filename][field] = status

    #     def analyze(self, data_txt):
    #         """Parse `fastqc_data.txt` file and store the results in a dictionary."""
    #         # Hypothetical implementation
    #         lines = data_txt.split('\n')
    #         for line in lines:
    #             parts = line.split('\t')
    #             if parts and parts[0] in data_fields:
    #                 self.data_dict[parts[0]] = parts[1:]

    def process_zipped_reports(self, directory):
        """Summarize FastQC reports in all zip files within a specified directory."""
        try:
            for root, _, files in os.walk(directory):
                for filename in files:
                    if filename.endswith("fastqc.zip"):
                        with zipfile.ZipFile(
                            os.path.join(root, filename), "r"
                        ) as zip_ref:
                            for file in zip_ref.namelist():
                                with zip_ref.open(file) as f:
                                    if "summary.txt" in file:
                                        self.summarize(f.readlines())
                                    elif "fastqc_data.txt" in file:
                                        print(f"Found data file: {file}")
        except Exception as e:
            print(f"Failed to summarize FastQC reports: {e}")

    def summarize_summaries(self, out_file):
        """Aggregate PASS/WARN/FAIL statuses for each sample and write to a file."""
        try:
            with open(out_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Sample"] + summary_fields)
                # sort the summaries by filename
                for filename, summary in sorted(self.summaries.items()):
                    writer.writerow(
                        [filename]
                        + [summary.get(field, "") for field in summary_fields]
                    )
        except Exception as e:
            print(f"Failed to write summary to file: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python qcsummary.py <directory>")
        sys.exit(1)

    directory = sys.argv[1]
    qc_manager = QCManager(directory)
    out_file = "qcsummary.csv"
    qc_manager.summarize_summaries(out_file)
