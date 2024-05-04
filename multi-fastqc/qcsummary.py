#! /usr/bin/env python3
import csv
import os
import zipfile


class FastQCModule:
    """A FastQC module found in the `fastqc_data.txt `file."""

    def __init__(self, name, status, headers, data):
        """Initialize a FastQC module with a name, headers, and data."""
        self.name = name
        self.status = status
        self.headers = headers
        self.module_data = data

    def __str__(self):
        """Return a formatted string representation of the FastQC module."""
        formatted_string = f"Module: {self.name} ({self.status})\n"
        formatted_string += f"Headers: {self.headers}\n"
        for data in self.module_data:
            formatted_string += f"{data}\n"
        return formatted_string


class QCReport:
    """An object indexable by its filename that contains the QC data modules."""

    def __init__(self, filename):
        """Initialize with a filename and an empty dictionary of modules."""
        self.filename = filename
        self.modules = {}
        self.performance_metrics = {}

    def parse_basic_metrics(self):
        """Parse basic metrics from the Basic Statistics module, if present."""
        basic_stats = self.modules.get("Basic Statistics")
        if basic_stats:
            for data in basic_stats.module_data:
                key, value = data[0], data[1]
                self.performance_metrics[data[0]] = data[1]

    def pm2csv(self, outfile):
        """Write the performance metrics to a CSV file."""
        try:
            with open(outfile, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(self.header)
                writer.writerow(
                    [self.performance_metrics.get(h, "N/A") for h in self.header]
                )
        except Exception as e:
            print(f"Failed to write performance metrics to CSV: {e}")

    def formatted_string(self):
        """Return a formatted string representation of the QC report."""
        formatted_string = f"Report for {self.filename}\n"
        for module_name, module in self.modules.items():
            formatted_string += str(module) + "\n\n"
        return formatted_string


class QCManager:
    def __init__(self, directory=None):
        self.summaries = {}
        self.reports = {}  # Using a dictionary to manage QCReport objects by filename

    def process_directory(self, directory):
        zipfiles = [f for f in os.listdir(directory) if f.endswith("_fastqc.zip")]
        for zip_path in zipfiles:
            report = QCReport(zip_path)
            self.reports[zip_path] = report
            try:
                with zipfile.ZipFile(os.path.join(directory, zip_path), "r") as z:
                    for fname in z.namelist():
                        with z.open(fname) as file:
                            if "summary.txt" in fname:
                                self.summarize(file, zip_path)
                            elif "fastqc_data.txt" in fname:
                                self.analyze(file, zip_path)
            except Exception as e:
                print(f"Failed to process file {zip_path}: {e}")

    def summarize(self, file, filename):
        lines = [line.decode().strip() for line in file]
        for line in lines:
            status, field = line.split("\t")[:2]
            self.summaries.setdefault(filename, {})[field] = status

    def analyze(self, file, filename):
        lines = [line.decode().strip() for line in file]
        if not lines[0].startswith("##FastQC"):
            raise ValueError("This is not a valid FastQC file.")
        current_module = None
        for line in lines[1:]:
            if line.startswith(">>"):
                if line.startswith(">>END_MODULE"):
                    current_module = None
                else:
                    m_name, m_status = line[2:].split("\t")
                    current_module = FastQCModule(m_name, m_status, [], [])
                    self.reports[filename].modules[m_name] = current_module
            elif line.startswith("#") and current_module:
                if line.startswith("#Total Deduplicated Percentage"):

                    self.reports[filename].performance_metrics[
                        "Total Deduplicated Percentage"
                    ] = line.split("\t")[1]
                else:
                    current_module.headers = line[1:].split("\t")
            elif current_module:
                current_module.module_data.append(line.split("\t"))
        self.reports[filename].parse_basic_metrics()

    def summary2csv(self, outfile):
        """Aggregate the summaries of all reports and write to a CSV file."""
        try:
            with open(outfile, "w", newline="") as f:
                writer = csv.writer(f)
                headers = ["Filename"] + sorted(
                    next(iter(self.summaries.values())).keys()
                )
                writer.writerow(headers)
                for filename, summary in sorted(self.summaries.items()):
                    row = [filename] + [
                        summary.get(field, "N/A") for field in headers[1:]
                    ]
                    writer.writerow(row)
        except Exception as e:
            print(f"Failed to write summary to CSV: {e}")

    def pm2csv(self, outfile):
        """Write the performance metrics of all reports to a CSV file."""
        try:
            with open(outfile, "w", newline="") as f:
                writer = csv.writer(f)
                first_report = next(iter(self.reports.values()))
                if first_report.performance_metrics:
                    headers = ["Sample ID"] + list(
                        first_report.performance_metrics.keys()
                    )
                    writer.writerow(headers)
                    for filename, report in sorted(self.reports.items()):
                        metrics = [
                            report.performance_metrics.get(h, "N/A")
                            for h in headers[1:]
                        ]
                        sample_id = filename.split("fastqc.zip")[0]
                        if sample_id.endswith("_001"):
                            sample_id = sample_id[:-4]
                        writer.writerow([sample_id] + metrics)
        except Exception as e:
            print(f"Failed to write performance metrics to CSV: {e}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python qcsummary.py <directory>")
        sys.exit(1)
    qc = QCManager(sys.argv[1])
    qc.process_directory(sys.argv[1])
    qc.summary2csv("summary.csv")
    qc.pm2csv("performance_metrics.csv")
