import os

# import time
import subprocess
import urllib.request
import threading
from pathlib import Path
from typing import Literal, Dict, Union, List
from .colorprinter import ColorPrinter as pr
from .progressbar import ProgressBar

# Constants
GENOMES_MIRROR = "ftp://ftp.ncbi.nlm.nih.gov"

# Define file types and suffixes
FILE_SPECS = {
    "bwa_index": "fna.bwa_index.tar.gz",
    "bowtie_index": "fna.bowtie_index.tar.gz",
    "samtools_index": "fna.fai",
    "fasta": "fna.gz",
    "hisat2_index": "fna.hisat2_index.tar.gz",
    "refseq_gff": "refseq_annotation.gff.gz",
    "refseq_gtf": "refseq_annotation.gtf.gz",
}

# Define species info
SPECIES_INFO: Dict[str, Dict[str, str]] = {
    "mouse": {
        "uri": "genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39",
        "ref": "GCA_000001635.9_GRCm39",
    },
    "human": {
        "uri": "genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38",
        "ref": "GCA_000001405.15_GRCh38",
    },
}

# Type aliases using Literal
Species = Literal[*SPECIES_INFO.keys()]
Files = Literal[*FILE_SPECS.keys()]


def _download_file(species: Species, file: Files, show_progress: bool) -> None:
    """Download a file for the given species and pipeline file with optional progress."""
    species_data = SPECIES_INFO[species]
    file_suffix = FILE_SPECS[file]
    file_url = f"{GENOMES_MIRROR}/{species_data['uri']}/seqs_for_alignment_pipelines.ucsc_ids/{species_data['ref']}_full_analysis_set.{file_suffix}"
    file_path = Path(f"~/cbmf/genomes/{species}/{file_suffix}").expanduser()

    pr.info(f"Downloading {species} {file} from {file_url}...")
    try:
        with urllib.request.urlopen(file_url) as response:
            total_size = int(response.info().get("Content-Length", -1))
            downloaded_size = 0

            with open(file_path, "wb") as out_file:
                progress_bar = (
                    ProgressBar(
                        total_size, prefix=f"{file:<15} ({total_size})", length=30
                    )
                    if show_progress and total_size > 0
                    else None
                )

                while True:
                    buffer = response.read(8192)
                    if not buffer:
                        break
                    downloaded_size += len(buffer)
                    out_file.write(buffer)

                    if progress_bar:
                        progress_bar.update(downloaded_size)

                # if progress_bar:
                #     progress_bar.finish()

        pr.success(f"Downloaded {species} {file} to {file_path}")

    except Exception as e:
        pr.error(f"Error downloading {file} for {species}: {str(e)}")
        raise


def download_helper(species: Species, files: Union[List[Files], str]) -> None:
    """Download files for a species, handling single, multiple, or all files."""
    if species not in SPECIES_INFO:
        raise ValueError(f"Invalid species: {species}")

    if isinstance(files, str):
        if files == "ALL":
            files = list(FILE_SPECS.keys())
        else:
            raise ValueError(
                "Invalid string input for files. Use 'ALL' to download all files or provide a list of files."
            )

    # Ensure 'files' is a list
    if not isinstance(files, list):
        raise ValueError("Invalid input for files parameter. Must be a list or 'ALL'.")

    # Validate files and download
    invalid_files = [file for file in files if file not in FILE_SPECS]
    if invalid_files:
        raise ValueError(f"Invalid file types: {', '.join(invalid_files)}")

    if len(files) == 1:
        # Handle single file download with progress
        _download_file(species, files[0], show_progress=True)
    else:
        # Handle multiple file downloads without progress
        threads = []
        for file in files:
            thread = threading.Thread(
                target=_download_file, args=(species, file, False)
            )
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

    pr.success(f"All specified downloads completed for {species}")


def decompress_files(path: str):
    """Try to decompress specific files in a directory."""
    files = os.listdir(path)

    for file in files:
        file_path = os.path.join(path, file)
        if file.endswith("index.tar.gz"):
            pr.info(f"Decompressing {file} with tar...")
            try:
                subprocess.run(["tar", "-xvzf", file_path, "-C", path], check=True)
                pr.success(f"Successfully extracted {file}")
            except subprocess.CalledProcessError as e:
                pr.error(f"Failed to extract {file}: {e}")
        elif file.endswith(".gz"):
            pr.info(f"Decompressing {file} with gunzip...")
            try:
                subprocess.run(["gunzip", "-v", file_path], check=True)
                pr.success(f"Successfully decompressed {file}")
            except subprocess.CalledProcessError as e:
                pr.error(f"Failed to decompress {file}: {e}")
        else:
            pr.info(f"Skipping {file}.")
