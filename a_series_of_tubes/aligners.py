import subprocess
import os

NUM_CORES = str(os.cpu_count())


def bowtie2_cmd(reference, r1, r2):
    return ["bowtie2", "-x", reference, "-1", r1, "-2", r2, "-p", NUM_CORES]


def hisat2_cmd(reference, r1, r2):
    return ["hisat2", "-x", reference, "-1", r1, "-2", r2, "-p", NUM_CORES, "--dta"]


def bwa_cmd(reference, r1, r2):
    return ["bwa", "mem", "-t", NUM_CORES, reference, r1, r2]


align_cmds = {
    "bowtie2": bowtie2_cmd,
    "hisat2": hisat2_cmd,
    "bwa": bwa_cmd,
}

index_cmds = {
    "bowtie2": lambda r: ["bowtie2-build", r, r],
    "hisat2": lambda r: ["hisat2-build", r, r],
    "bwa": lambda r: ["bwa", "index", r],
}


def align_reads(aligner, r1_fastq, r2_fastq, reference):
    cmd = align_cmds.get(aligner)
    if not cmd:
        raise ValueError(f"Invalid aligner: {aligner}")
    return subprocess.run(
        cmd(reference, r1_fastq, r2_fastq), capture_output=True, text=True, check=True
    ).stdout


def build_index(aligner, reference):
    cmd = index_cmds.get(aligner)
    if not cmd:
        raise ValueError(f"Unsupported aligner: {aligner}")
    subprocess.run(cmd(reference), check=True)
