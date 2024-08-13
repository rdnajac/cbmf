import os
import subprocess
from typing import Literal

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
) -> None:
    """
    Download a specific genome file for a given species.

    This function downloads a specified genome file for either mouse or human
    from the NCBI FTP server. It handles both direct downloads and
    download-and-extract operations for index files.

    :param species: The species for which to download the genome file.
    :type species: Literal['mouse', 'human']
    :param file: The type of file to download.
    :type file: Literal['bwa_index', 'samtools_index', 'fasta', 'hisat2_index', 'refseq_gff', 'refseq_gtf']
    :return: None
    :raises ValueError: If an invalid species or file type is provided.
    :raises subprocess.CalledProcessError: If the download or extraction process fails.

    :Example:

    >>> download_genome_file('mouse', 'bwa_index')
    Downloading bwa index for mouse...
    Download complete.
    """
    if species not in SPECIES_INFO:
        raise ValueError(f"Invalid species: {species}")

    if file not in PIPELINE_FILES:
        raise ValueError(f"Invalid file type: {file}")

    species_data = SPECIES_INFO[species]
    file_suffix = PIPELINE_FILES[file]
    file_url = f"{GENOMES_MIRROR}/{species_data['uri']}/seqs_for_alignment_pipelines.ucsc_ids/{species_data['ref']}_full_analysis_set.{file_suffix}"
    genomes_dir = os.path.expanduser(f"~/cbmf/genomes/{species}")
    os.makedirs(genomes_dir, exist_ok=True)

    try:
        subprocess.run(["curl", "-O", file_url], cwd=genomes_dir, check=True)
        # TODO decompress files

    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(f"Error downloading {file}: {str(e)}")


def download_all_files_for_species(species: Literal["mouse", "human"]) -> None:
    """
    Download all genome files for a given species.

    This function downloads all specified genome files for either mouse or human
    from the NCBI FTP server.

    :param species: The species for which to download all genome files.
    :type species: Literal['mouse', 'human']
    :return: None
    :raises ValueError: If an invalid species is provided.
    """
    if species not in SPECIES_INFO:
        raise ValueError(f"Invalid species: {species}")

    for file in PIPELINE_FILES.keys():
        try:
            print(f"Downloading {file} for {species}...")
            download_genome_file(species, file)
            print(f"Successfully downloaded {file} for {species}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {file} for {species}: {str(e)}")
        except Exception as e:
            print(f"Unexpected error downloading {file} for {species}: {str(e)}")
