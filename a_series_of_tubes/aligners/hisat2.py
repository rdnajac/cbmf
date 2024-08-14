import subprocess


def align_hisat2(r1_fastq, r2_fastq, reference):
    """
    Align paired-end reads using HISAT2.

    :param r1_fastq: Path to the first read file (R1)
    :param r2_fastq: Path to the second read file (R2)
    :param reference: Path to the reference genome
    :return: stdout buffer from the HISAT2 command
    """
    print(f"__{align_hisat2.__name__}__")
    return None
    cmd = [
        "hisat2",
        "-p",
        str(subprocess.os.cpu_count()),
        "--mm",
        "-1",
        r1_fastq,
        "-2",
        r2_fastq,
        "-x",
        reference,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout


# TODO: add function to build index
# hisat2-build genome.fa genome
