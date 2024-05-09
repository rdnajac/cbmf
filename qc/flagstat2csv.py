#! /usr/bin/env python3
import re
import sys
import os
import glob
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed


class Flagstat:
    def __init__(self, flagstat_output):
        self.data = flagstat_output
        self.stats = self.parse_stats()

    def parse_stats(self):
        patterns = {
            "total_reads": r"(\d+) \+ (\d+) in total",
            "mapped": r"(\d+) \+ (\d+) mapped",
            "duplicates": r"(\d+) \+ (\d+) duplicates",
        }
        results = {}
        for key, pattern in patterns.items():
            match = re.search(pattern, self.data)
            if match:
                results[key] = {
                    "qc_passed": int(match.group(1)),
                    "qc_failed": int(match.group(2)),
                }
        return results

    def calculate_additional_stats(self, total_reads):
        mapped_reads = self.stats.get("mapped", {}).get("qc_passed", 0)
        duplicate_reads = self.stats.get("duplicates", {}).get("qc_passed", 0)
        unique_reads = mapped_reads - duplicate_reads
        self.stats["% Mapped Reads"] = (
            (mapped_reads / total_reads) * 100 if total_reads else 0
        )
        self.stats["% Duplicate Reads"] = (
            (duplicate_reads / total_reads) * 100 if total_reads else 0
        )
        self.stats["Unique Reads"] = unique_reads
        self.stats["% Unique Reads"] = (
            (unique_reads / total_reads) * 100 if total_reads else 0
        )
        self.stats["% Unique of Mapped"] = (
            (unique_reads / mapped_reads) * 100 if mapped_reads else 0
        )


def process_bam_file(file):
    flagstat_output = os.popen(f"samtools flagstat {file}").read()
    flagstat = Flagstat(flagstat_output)
    total_reads = int(os.popen(f"samtools view -c {file}").read().strip())
    flagstat.calculate_additional_stats(total_reads)
    return os.path.basename(file), flagstat.stats


if __name__ == "__main__":
    directory = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    bam_files = glob.glob(os.path.join(directory, "*.bam"))
    results_dict = {}

    with ThreadPoolExecutor(max_workers=64) as executor:
        futures = {executor.submit(process_bam_file, file): file for file in bam_files}
        for future in as_completed(futures):
            filename, stats = future.result()
            results_dict[filename] = stats

    headers = [
        "Sample ID",
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

    with open("flagstat.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for filename, stats in results_dict.items():
            total_reads = stats.get("total_reads", {}).get("qc_passed", 0) + stats.get(
                "total_reads", {}
            ).get("qc_failed", 0)
            mapped_reads = stats.get("mapped", {}).get("qc_passed", 0)
            duplicate_reads = stats.get("duplicates", {}).get("qc_passed", 0)
            unique_reads = stats.get("Unique Reads", 0)
            row = [
                filename.split(".")[0],
                total_reads,
                stats.get("total_reads", {}).get("qc_passed", 0),
                (
                    100 * stats.get("total_reads", {}).get("qc_passed", 0) / total_reads
                    if total_reads
                    else 0
                ),
                mapped_reads,
                stats.get("% Mapped Reads", 0),
                duplicate_reads,
                stats.get("% Duplicate Reads", 0),
                unique_reads,
                stats.get("% Unique Reads", 0),
                stats.get("% Unique of Mapped", 0),
            ]
            writer.writerow(row)
