#!/usr/bin/env python3
"""
MD5Handler: A script to compute and manage MD5 hashes for files in a directory.

This script uses the system's `md5sum` command to generate and verify MD5 hashes.
"""

import os
import subprocess
import sys


class MD5Handler:
    """
    A class to handle MD5 hash computation and validation for files in a directory.
    """

    def __init__(self, directory, extension=".fastq.gz", output_file="md5sums.txt"):
        """
        Initialize the MD5Handler.

        Args:
            directory (str): Path to the directory to process.
            extension (str): File extension to filter files. Defaults to ".fastq.gz".
            output_file (str): Name of the file to write consolidated MD5 hashes. Defaults to "md5sums.txt".
        """
        self.directory = directory
        self.extension = extension
        self.output_file = output_file
        self.md5s = {}
        self.files = []

        self._validate_directory()
        self._get_files()
        self._compute_md5s()

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
            print(
                f"No files with extension '{self.extension}' found in {self.directory}."
            )

    def _compute_md5s(self):
        """Compute MD5 hashes for files and store them in a dictionary."""
        for file in self.files:
            file_path = os.path.join(self.directory, file)
            try:
                result = subprocess.run(
                    ["md5sum", file_path], capture_output=True, text=True, check=True
                )
                self.md5s[file] = result.stdout.split()[0]
            except subprocess.CalledProcessError as e:
                print(f"Error computing MD5 for {file}: {e}")
            except Exception as e:
                print(f"Unexpected error with {file}: {e}")

    def check_md5s(self):
        """
        Recompute and verify MD5 hashes for files.

        Prints a message if there are mismatches.
        """
        print("Checking MD5 hashes...")
        for file in self.files:
            file_path = os.path.join(self.directory, file)
            try:
                result = subprocess.run(
                    ["md5sum", file_path], capture_output=True, text=True, check=True
                )
                computed_md5 = result.stdout.split()[0]
                if computed_md5 != self.md5s[file]:
                    print(f"MD5 mismatch for {file}")
            except subprocess.CalledProcessError as e:
                print(f"Error verifying MD5 for {file}: {e}")
            except Exception as e:
                print(f"Unexpected error with {file}: {e}")

    def consolidate_md5s(self):
        """
        Write computed MD5 hashes to an output file.

        Each line in the file will contain the MD5 hash and the corresponding filename.
        """
        output_path = os.path.join(self.directory, self.output_file)
        try:
            with open(output_path, "w") as f:
                for file, md5 in self.md5s.items():
                    f.write(f"{md5}  {file}\n")
            print(f"MD5 hashes written to {output_path}.")
        except IOError as e:
            print(f"Error writing to {output_path}: {e}")


# Entry point for the script
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 md5handler.py <directory>")
        sys.exit(1)

    directory = sys.argv[1]

    try:
        handler = MD5Handler(directory)
        handler.check_md5s()
        handler.consolidate_md5s()
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
