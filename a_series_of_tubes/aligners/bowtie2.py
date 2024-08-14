import subprocess


def align_bowtie2(r1_fastq, r2_fastq, reference):
    """
    Align paired-end reads using Bowtie2.

    :param r1_fastq: Path to the first read file (R1)
    :param r2_fastq: Path to the second read file (R2)
    :param reference: Path to the reference genome
    :return: stdout buffer from the Bowtie2 command
    """
    print(f"__{align_bowtie2.__name__}__")
    return None
    cmd = [
        "bowtie2",
        "-p",
        str(subprocess.os.cpu_count()),
        "--mm",
        "-1",
        r1_fastq,
        "-2",
        r2_fastq,
        "-x",
        reference,
        "--dta",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout


# TODO: add function to build index
# bowtie2-build example/reference/lambda_virus.fa example/index/lambda_virus
# bowtie2-build --large-index example/reference/lambda_virus.fa example/index/lambda_virus
