import subprocess
from ..utils.logger import logger


def align_PE_reads(aligner, reference, r1fastq, r2fastq, threads):
    """Align paired-end reads to a reference genome."""
    # TODO
    pass


def align_reads(aligner, input_dir, PE=True, reference=None, species=None, threads=1):
    """
    Align reads to reference genome

    :param aligner: str, aligner to use (bwa, bowtie2, hisat2)
    :param input_dir: Path, input directory containing FASTQ files
    :param PE: bool, whether reads are paired-end
    :param reference: Path, path to reference genome index (optional if species is provided)
    :param species: str, species name (optional if reference is provided)
    :param threads: int, number of threads to use
    """
    if not PE:
        raise ValueError("Only paired-end reads are supported")

    logger.info(f"Starting read alignment with {aligner}")

    if reference is None and species is None:
        raise ValueError("Either reference or species must be provided")
    elif reference is None:
        from .genomemanager import GenomeManager

        reference = GenomeManager.get_genome_index(species, aligner)
    elif reference.exists():
        logger.info(f"Using reference genome index: {reference}")
    else:
        raise FileNotFoundError(f"Reference genome index not found: {reference}")

    align_PE_reads_from_dir(aligner, input_dir, reference, threads)


def align_PE_reads_from_dir(aligner, input_dir, reference, threads):

    # Find all *_R1.fastq.gz files and their matching *_R2.fastq.gz files
    r1fastqs = sorted(input_dir.glob("*_R1.fastq.gz"))
    paired_files = []

    for r1fastq in r1fastqs:
        r2fastq = r1fastq.with_name(
            r1fastq.name.replace("_R1.fastq.gz", "_R2.fastq.gz")
        )
        if r2fastq.exists():
            paired_files.append((r1fastq, r2fastq))
        else:
            raise FileNotFoundError(f"No matching R2 file found for {r1fastq.name}")

    if not paired_files:
        raise FileNotFoundError(
            f"No paired FASTQ files found in input directory: {input_dir}"
        )

    for r1fastq, r2fastq in paired_files:
        if aligner == "bwa":
            raise NotImplementedError("BWA alignment not implemented yet")
        elif aligner == "bowtie2":
            _bowtie2(r1fastq, r2fastq, reference, threads)
        elif aligner == "hisat2":
            _hisat2(r1fastq, r2fastq, reference, threads)


def _bwa(reference, r1fastq, r2fastq, output_file, threads):
    cmd = [
        "bwa",
        "mem",
        "-t",
        threads,
        str(reference),
        str(r1fastq),
        str(r2fastq),
    ]
    _run_command(cmd)


def _bowtie2(reference, r1fastq, r2fastq, threads):
    cmd = [
        "bowtie2",
        "-x",
        str(reference),
        "-1",
        str(r1fastq),
        "-2",
        str(r2fastq),
        "-p",
        threads,
    ]
    _run_command(cmd)


def _hisat2(reference, r1fastq, r2fastq, threads):
    cmd = [
        "hisat2",
        "-x",
        str(reference),
        "-1",
        str(r1fastq),
        "-2",
        str(r2fastq),
        "-p",
        threads,
    ]
    _run_command(cmd)


def _run_command(cmd):
    try:
        subprocess.run(" ".join(cmd), shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Error running command: {' '.join(cmd)}") from e
