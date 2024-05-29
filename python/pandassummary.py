import os
import re
import pandas as pd
import zipfile
import sys

csv_headers = [
    "Sample ID",
    "Yield (Mbases)",
    "Mean Quality Score",
    "% Bases >= 30",
    "[paired-end] Reads",
    "Mapped Reads",
    "% Mapped Reads",
    "Duplicate Reads",
    "% Duplicate Reads (of all reads)",
    "Unique Reads",
    "% Unique Reads (of all reads)",
    "% Unique (of mapped)",
]




class FastQCAnalyzerWithPandas:
    def __init__(self, directory):
        self.directory = directory
        self.summary_df = pd.DataFrame(columns=csv_headers)

    def parse_file_content(self, content, file_name):
        nseq = re.search(r"Total Sequences\t(\d+)", content).group(1)
        dup = re.search(r"Total Deduplicated Percentage\t(\d+\.\d)", content).group(1)
        kmer = re.search(r"Kmer Content\t(\w+)", content).group(1)
        adapter = re.search(r"Adapter Content\t(\w+)", content).group(1)
        gc_end = self.calculate_gc_skew(content)
        sample, rep, ge, lane, direction = self.get_lib_info(file_name)
        self.summary_df = self.summary_df.append(
            {
                "Sample ID": sample,
                "Replicate": rep,
                "Group": ge,
                "Lane": lane,
                "Direction": direction,
                "Total Sequences": nseq,
                "Deduplicated Percentage": dup,
                "GC Skew": gc_end,
                "Kmer Content": kmer,
                "Adapter Content": adapter,
            },
            ignore_index=True,
        )

    def calculate_gc_skew(self, content):
        gc_values = re.search(r"\n145-149\s([0-9\.\s]+)\n", content).group(1).split()
        gc, at = float(gc_values[0]), float(gc_values[3])
        return str(round(gc / (gc + at), 2))

    def get_lib_info(self, file_name):
        match = re.match(
            r"Anc-(?P<sample>\d+)-(?P<rep>[AB])_(?P<GE>S\d+)_(?P<lane>L\d\d\d)_(?P<direction>R\d)_fastqc",
            file_name,
        )
        return match.groups()

    def process_single_file(self, file_name):
        if file_name.endswith(".txt"):
            with open(os.path.join(self.directory, file_name), "r") as file:
                content = file.read()
                self.parse_file_content(content, file_name)

    def process_files(self):
        for file_name in os.listdir(self.directory):
            self.process_single_file(file_name)
        self.summary_df = self.summary_df.reindex(columns=csv_headers)

    def print_single_summary(self, sample_id):
        print(self.summary_df[self.summary_df["Sample ID"] == sample_id])

    def print_summary(self):
        print(self.summary_df)


def parse_stats(lines, summary, detailed_summary):
    regex = re.compile(r"(\w+)\s+(.+?)\s+(\S+\.fastq\.gz)")
    for line in lines:
        match = regex.match(line.strip())
        if match:
            status, test_name, file_name = match.groups()
            summary[status] = summary.get(status, 0) + 1
            if status == "FAIL" or status == "WARN":
                detailed_summary.setdefault(file_name, []).append((status, test_name))


def multireport(directory):
    os.chdir(directory)
    report_filename = "multiqc_report.txt"
    files = [f for f in os.listdir(".") if f.endswith("fastqc.zip")]
    summary = {}
    detailed_summary = {}
    if not files:
        print("No fastqc.zip files found", file=sys.stderr)
        return
    with open(report_filename, "w") as report_file:
        for file in files:
            with zipfile.ZipFile(file, "r") as zip_ref:
                summary_files = [f for f in zip_ref.namelist() if "summary.txt" in f]
                for summary_file in summary_files:
                    with zip_ref.open(summary_file) as f:
                        lines = f.read().decode("utf-8").splitlines()
                        report_file.writelines(line + "\n" for line in lines)
                        report_file.write("----------------------------------------\n")
                        parse_stats(lines, summary, detailed_summary)
    print_summary(summary, detailed_summary)


def print_summary(summary, detailed_summary):
    statuses = ["PASS", "FAIL", "WARN"]
    print("\nFastQC Summary Statistics:", file=sys.stderr)
    for status in statuses:
        count = summary.get(status, 0)
        print(f"{status}: {count}", file=sys.stderr)
    print("\nDetailed Fail and Warn Statistics:", file=sys.stderr)
    for file_name, issues in detailed_summary.items():
        for issue in issues:
            print(f"{file_name}: {issue[0]} - {issue[1]}", file=sys.stderr)


# Example usage
analyzer = FastQCAnalyzerWithPandas("./")
analyzer.process_files()
analyzer.print_summary()
multireport("your_directory_path")
