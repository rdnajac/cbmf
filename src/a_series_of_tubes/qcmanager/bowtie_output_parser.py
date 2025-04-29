#!/usr/bin/env python3

import csv
import os
import re


class BowtieOutputDictionary:
    def __init__(self):
        self.samples = []
        self.bowtie_output = {}

    def add_sample(self, sample_name, sample_output):
        self.samples.append(sample_name)
        self.bowtie_output[sample_name] = sample_output

    def read_logs(self, folder):
        for file in os.listdir(folder):
            if file.endswith(".log"):
                sample_name = file[:-4]
                with open(os.path.join(folder, file), "r") as f:
                    self.add_sample(sample_name, f.read())

    def parse_output(self):
        parsed_data = []
        for sample in self.samples:
            output = self.bowtie_output[sample]
            total_reads = int(re.search(r"(\d+) reads;", output).group(1))
            aligned_once = int(
                re.search(
                    r"(\d+) \(\d+\.\d+%\) aligned concordantly exactly 1 time", output
                ).group(1)
            )
            aligned_more_than_once = int(
                re.search(
                    r"(\d+) \(\d+\.\d+%\) aligned concordantly >1 times", output
                ).group(1)
            )
            overall_alignment_rate = float(
                re.search(r"(\d+\.\d+)% overall alignment rate", output).group(1)
            )

            aligned_reads = aligned_once + aligned_more_than_once

            parsed_data.append(
                {
                    "Sample": sample,
                    "Total Reads": total_reads,
                    "Aligned Reads": aligned_reads,
                    "Overall Alignment Rate (%)": overall_alignment_rate,
                }
            )
        return parsed_data

    def write_to_csv(self, output_csv, parsed_data):
        fieldnames = [
            "Sample",
            "Total Reads",
            "Aligned Reads",
            "Overall Alignment Rate (%)",
        ]

        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for data in parsed_data:
                writer.writerow(data)


def combine_run_csvs(run_folders, output_csv):
    all_data = []
    for run_folder in run_folders:
        bowtie_data = BowtieOutputDictionary()
        bowtie_data.reata.parse_output()

    bowtie_data.write_to_csv(output_csv, all_data)


# Example usage
base_folder = "../data/kalay/qc/"
run_folders = [os.path.join(base_folder, f"run{i}/bowtie2log") for i in range(1, 6)]
combined_output_csv = "~/cbmf/data/kalay/qc/bowtie2_summary_all_runs.csv"
combine_run_csvs(run_folders, combined_output_csv)
