import os
import sys
import time
import urllib.request
import threading
import subprocess
from typing import Literal, Dict
from .colorprinter import ColorPrinter as pr
from .progressbar import ProgressBar

# Constants
GENOMES_MIRROR = "ftp://ftp.ncbi.nlm.nih.gov"

# Configuration Data
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

PIPELINE_FILES: Dict[str, str] = {
    "bwa_index": "fna.bwa_index.tar.gz",
    "bowtie_index": "fna.bowtie_index.tar.gz",
    "samtools_index": "fna.fai",
    "fasta": "fna.gz",
    "hisat2_index": "fna.hisat2_index.tar.gz",
    "refseq_gff": "refseq_annotation.gff.gz",
    "refseq_gtf": "refseq_annotation.gtf.gz",
}

# Type aliases
Species = Literal["mouse", "human"]
PipelineFile = Literal[
    "bwa_index", "samtools_index", "fasta", "hisat2_index", "refseq_gff", "refseq_gtf"
]


def download_genome_file(
    species: Species, file: PipelineFile, show_progress: bool = True
) -> None:
    if species not in SPECIES_INFO:
        raise ValueError(f"Invalid species: {species}")
    if file not in PIPELINE_FILES:
        raise ValueError(f"Invalid file type: {file}")

    species_data = SPECIES_INFO[species]
    file_suffix = PIPELINE_FILES[file]
    file_url = f"{GENOMES_MIRROR}/{species_data['uri']}/seqs_for_alignment_pipelines.ucsc_ids/{species_data['ref']}_full_analysis_set.{file_suffix}"
    genomes_dir = os.path.expanduser(f"~/cbmf/genomes/{species}")
    os.makedirs(genomes_dir, exist_ok=True)

    pr.info(f"Downloading {file} for {species} from {file_url}...")

    try:
        with urllib.request.urlopen(file_url) as response:
            total_size = int(response.info().get("Content-Length", -1))
            downloaded_size = 0
            file_path = os.path.join(genomes_dir, os.path.basename(file_url))

            last_update_time = 0

        with open(file_path, "wb") as out_file:
            progress_bar = ProgressBar(
                total_size, prefix=f"{file:<15} ({total_size})", length=30
            )
            while True:
                buffer = response.read(8192)
                if not buffer:
                    break
                downloaded_size += len(buffer)
                out_file.write(buffer)
                if show_progress and total_size > 0:
                    current_time = time.time()
                    if current_time - last_update_time >= 2:
                        progress_bar.update(downloaded_size)
                        last_update_time = current_time

                progress_bar.finish()
                pr.success(f"Successfully downloaded {file} for {species}")

    except Exception as e:
        pr.error(f"Error downloading {file} for {species}: {str(e)}")
        raise


def download_all_files_for_species(species: Species) -> None:
    if species not in SPECIES_INFO:
        raise ValueError(f"Invalid species: {species}")

    threads = []
    for file in PIPELINE_FILES.keys():
        thread = threading.Thread(
            target=download_genome_file, args=(species, file, False)
        )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    pr.success(f"All downloads completed for {species}")


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
