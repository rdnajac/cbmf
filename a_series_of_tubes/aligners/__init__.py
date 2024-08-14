from .bowtie2 import align_bowtie2
from .hisat2 import align_hisat2
from .bwa import align_bwa


def align_reads(aligner, r1_fastq, r2_fastq, reference):
    """
    Generic function to call specific aligner functions.

    :param aligner: Name of the aligner to use ('bowtie2' or 'hisat2')
    :param r1_fastq: Path to the first read file (R1)
    :param r2_fastq: Path to the second read file (R2)
    :param reference: Path to the reference genome
    :return: stdout buffer from the aligner command
    """
    print(f'align reads with args: {aligner}, {r1_fastq}, {r2_fastq}, {reference}')
    if aligner == "bowtie2":
        return align_bowtie2(r1_fastq, r2_fastq, reference)
    elif aligner == "hisat2":
        return align_hisat2(r1_fastq, r2_fastq, reference)
    elif aligner == "bwa":
        return align_bwa(r1_fastq, r2_fastq, reference)
    elif aligner == "STAR":
        print(f"{aligner} not implemented yet")
        return None
    else:
        raise ValueError(f"Invalid aligner: {aligner}")
