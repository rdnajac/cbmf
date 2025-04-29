#!/usr/bin/env python3
"""
MD5Handler: A script to validate and manage MD5 hashes for files in a directory.

This script checks for individual MD5 files corresponding to each data file
and validates their integrity using multithreading.
"""

import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


class MD5Handler:
    """
    A class to handle MD5 validation and management for files in a directory.
    """

    def __init__(
        self, directory, extension=".fastq.gz", output_file="md5sums.txt", threads=4
    ):
        """
        Initialize the MD5Handler.

        Args:
            directory (str): Path to the directory to process.
            extension (str): File extension to filter files. Defaults to ".fastq.gz".
            output_file (str): Name of the file to write consolidated MD5 hashes. Defaults to "md5sums.txt".
            threads (int): Number of threads to use for parallel processing. Defaults to 4.
        """
        self.directory = directory
        self.extension = extension
        self.output_file = output_file
        self.threads = threads
        self.files = []
        self.md5_files = []

        self._validate_directory()
        self._get_files()
        self._get_md5_files()

    def _validate_directory(self):
        """Validate that the given directory exists and is accessible."""
        if not os.path.isdir(self.directory):
            raise NotADirectoryError(
                f"Error: {self.directory} is not a valid directory."
            )

    def _get_files(self):
        """Retrieve files in the directory matching the specified extension."""
        self.files = [
            file
            for file in os.listdir(self.directory)
            if file.endswith(self.extension)
            and os.path.isfile(os.path.join(self.directory, file))
        ]
        if not self.files:
            raise FileNotFoundError(
                f"No files with extension '{self.extension}' found in {self.directory}."
            )

    def _get_md5_files(self):
        """Retrieve MD5 files corresponding to the data files."""
        self.md5_files = [
            f"{file}.md5"
            for file in self.files
            if os.path.isfile(os.path.join(self.directory, f"{file}.md5"))
        ]
        missing_md5 = set(self.files) - set(file[:-4] for file in self.md5_files)
        if missing_md5:
            raise FileNotFoundError(f"Missing MD5 files for: {', '.join(missing_md5)}")

    def _validate_md5(self, file):
        """
        Validate the MD5 hash of a single file.

        Args:
            file (str): The file to validate.

        Raises:
            ValueError: If the MD5 hash does not match.
        """
        file_path = os.path.join(self.directory, file)
        md5_file_path = os.path.join(self.directory, f"{file}.md5")
        try:
            with open(md5_file_path, "r") as md5_file:
                expected_md5 = md5_file.read().strip().split()[0]
            result = subprocess.run(
                ["md5sum", file_path], capture_output=True, text=True, check=True
            )
            computed_md5 = result.stdout.split()[0]
            if computed_md5 != expected_md5:
                raise ValueError(
                    f"MD5 mismatch for {file}: expected {expected_md5}, got {computed_md5}"
                )
            return (file, expected_md5)
        except Exception as e:
            return (file, str(e))

    def validate_md5s(self):
        """Validate MD5 hashes for all files in parallel."""
        print(f"Validating MD5 hashes with {self.threads} threads...")
        results = []
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            future_to_file = {
                executor.submit(self._validate_md5, file): file for file in self.files
            }
            for future in as_completed(future_to_file):
                results.append(future.result())

        # Check results for errors
        errors = [
            result
            for result in results
            if isinstance(result[1], str) and "MD5 mismatch" in result[1]
        ]
        if errors:
            for error in errors:
                print(error[1])
            sys.exit(1)
        print("All MD5 hashes validated successfully.")

    def consolidate_md5s(self):
        """
        Write validated MD5 hashes to an output file.

        Each line in the file will contain the MD5 hash and the corresponding filename.
        """
        output_path = os.path.join(self.directory, self.output_file)
        try:
            with open(output_path, "w") as f:
                for file in self.files:
                    md5_file_path = os.path.join(self.directory, f"{file}.md5")
                    with open(md5_file_path, "r") as md5_file:
                        md5 = md5_file.read().strip().split()[0]
                    f.write(f"{md5}  {file}\n")
            print(f"Validated MD5 hashes consolidated into {output_path}.")
        except IOError as e:
            print(f"Error writing to {output_path}: {e}")


# Entry point for the script
if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python3 md5handler.py <directory> [threads]")
        sys.exit(1)

    directory = sys.argv[1]
    threads = int(sys.argv[2]) if len(sys.argv) == 3 else 4

    try:
        handler = MD5Handler(directory, threads=threads)
        handler.validate_md5s()
        handler.consolidate_md5s()
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
