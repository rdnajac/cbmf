from pathlib import Path

# FTP server for genome downloads
GENOMES_MIRROR = "ftp://ftp.ncbi.nlm.nih.gov"

# Reference genome information
REFERENCE = {
    "mouse": {
        "uri": "genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39",
        "ref": "GCA_000001635.9_GRCm39",
    },
    "human": {
        "uri": "genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38",
        "ref": "GCA_000001405.15_GRCh38",
    },
}

# File types available for download
FILES = {
    "bwa_index": "fna.bwa_index.tar.gz",
    "bowtie_index": "fna.bowtie_index.tar.gz",
    "samtools_index": "fna.fai",
    "fasta": "fna.gz",
    "hisat2_index": "fna.hisat2_index.tar.gz",
    "refseq_gff": "refseq_annotation.gff.gz",
    "refseq_gtf": "refseq_annotation.gtf.gz",
}

# Path to the scripts directory
SCRIPTS_DIR = Path.home() / "cbmf" / "scripts"

# Add any other global constants here
