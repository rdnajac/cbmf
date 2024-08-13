# utils/genomemanager.py

import os
import urllib.request
import threading
from typing import Literal
from .colorprinter import ColorPrinter as pr

GENOMES_MIRROR = "ftp://ftp.ncbi.nlm.nih.gov"
SPECIES_INFO = {
    "mouse": {
        "uri": "genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39",
        "ref": "GCA_000001635.9_GRCm39",
    },
    "human": {
        "uri": "genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38",
        "ref": "GCA_000001405.15_GRCh38",
    },
}
PIPELINE_FILES = {
    "bwa_index": "fna.bwa_index.tar.gz",
    "samtools_index": "fna.fai",
    "fasta": "fna.gz",
    "hisat2_index": "fna.hisat2_index.tar.gz",
    "refseq_gff": "refseq_annotation.gff.gz",
    "refseq_gtf": "refseq_annotation.gtf.gz",
}

def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█', print_end="\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
    if iteration == total: 
        print()

def download_genome_file(
    species: Literal["mouse", "human"],
    file: Literal[
        "bwa_index",
        "samtools_index",
        "fasta",
        "hisat2_index",
        "refseq_gff",
        "refseq_gtf",
    ],
    show_progress: bool = True
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
    
    pr.info(f"Downloading {file} for {species}...")
    
    try:
        with urllib.request.urlopen(file_url) as response:
            total_size = int(response.info().get('Content-Length', -1))
            downloaded_size = 0
            
            with open(os.path.join(genomes_dir, os.path.basename(file_url)), 'wb') as out_file:
                while True:
                    buffer = response.read(8192)
                    if not buffer:
                        break
                    downloaded_size += len(buffer)
                    out_file.write(buffer)
                    if show_progress:
                        print_progress_bar(downloaded_size, total_size, prefix=f'{file:<15}', length=30)
        
        pr.success(f"Successfully downloaded {file} for {species}")
        # TODO: decompress files
    
    except Exception as e:
        pr.error(f"Error downloading {file} for {species}: {str(e)}")
        raise

def download_all_files_for_species(species: Literal["mouse", "human"]) -> None:
    if species not in SPECIES_INFO:
        raise ValueError(f"Invalid species: {species}")
    
    threads = []
    lock = threading.Lock()

    for file in PIPELINE_FILES.keys():
        thread = threading.Thread(target=download_genome_file, args=(species, file, False))
        threads.append(thread)
        thread.start()
    
    for thread in threads:
        thread.join()
    
    pr.success(f"All downloads completed for {species}")

def main():
    download_genome_file("mouse", "bwa_index")

if __name__ == "__main__":
    main()
