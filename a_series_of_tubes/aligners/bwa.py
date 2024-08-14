import os
import subprocess


def align_bwa(r1_fastq, r2_fastq, reference):
    """
    Align paired-end reads using BWA-MEM.

    :param r1_fastq: Path to the first read file (R1)
    :param r2_fastq: Path to the second read file (R2)
    :param reference: Path to the reference genome
    :return: stdout buffer from the BWA command
    """
    print(f"__{align_bwa.__name__}__")
    return None
    cmd = [
        "bwa",
        "mem",
        "-t",
        str(os.cpu_count()),  # Use all available CPU cores
        "-M",  # Mark shorter split hits as secondary
        reference,
        r1_fastq,
        r2_fastq,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout


def index_reference(reference):
    """
    Index the reference genome using BWA.

    :param reference: Path to the reference genome
    """
    print(f"__{index_reference.__name__}__")
    return None
    cmd = ["bwa", "index", reference]
    subprocess.run(cmd, check=True)
