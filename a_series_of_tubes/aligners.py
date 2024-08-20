import subprocess
import os

NUM_CORES = str(os.cpu_count())

align_cmds = {
    "bowtie2": lambda r, r1, r2: ["bowtie2", "-x", r, "-1", r1, "-2", r2, "-p", NUM_CORES],
    "hisat2": lambda r, r1, r2: ["hisat2", "-x", r, "-1", r1, "-2", r2, "-p", NUM_CORES, "--dta"],
    "bwa": lambda r, r1, r2: ["bwa", "mem", "-t", NUM_CORES, r, r1, r2],
    "star": lambda r, r1, r2: ["STAR", "--genomeDir", r, "--readFilesIn", r1, r2, "--runThreadN", NUM_CORES]
}

index_cmds = {
    "bowtie2": lambda r: ["bowtie2-build", r, r],
    "hisat2": lambda r: ["hisat2-build", r, r],
    "bwa": lambda r: ["bwa", "index", r],
    "star": lambda r: ["STAR", "--runMode", "genomeGenerate", "--genomeDir", os.path.dirname(r), "--genomeFastaFiles", r]
}

def align_reads(aligner, r1_fastq, r2_fastq, reference):
    cmd = align_cmds.get(aligner)
    if not cmd:
        raise ValueError(f"Invalid aligner: {aligner}")
    return subprocess.run(cmd(reference, r1_fastq, r2_fastq), capture_output=True, text=True, check=True).stdout

def build_index(aligner, reference):
    cmd = index_cmds.get(aligner)
    if not cmd:
        raise ValueError(f"Unsupported aligner: {aligner}")
    subprocess.run(cmd(reference), check=True)

def main(aligner, r1_fastq, r2_fastq, reference):
    build_index(aligner, reference)
    alignment_output = align_reads(aligner, r1_fastq, r2_fastq, reference)
    print(alignment_output)

if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
