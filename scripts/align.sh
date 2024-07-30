#!/bin/bash
#
## Align a pair of fastq files using HISAT2 or Bowtie2.
set -euo pipefail
set -x

MOUSEREF="/home/ubuntu/genomes/mouse/GCA_000001635.9_GRCm39_full_analysis_set.fna.hisat2_index"
NUMTHREADS="$(nproc)"

# expects to read in from stdin
# finally pipe to stdout and leave it to the caller to redirect to a file
set -x
export NUMTHREADS=16
samtools_to_bam()
{
	samtools view -@ "$NUMTHREADS" $1 -b \
		| samtools sort -n -@ "$NUMTHREADS" - \
		| samtools fixmate -@ "$NUMTHREADS" -m - - \
		| samtools sort -@ "$NUMTHREADS" - \
		| samtools view -b -@ "$NUMTHREADS" -
}

for sam_file in ./*.sam; do
	samtools_to_bam "$sam_file" > "${sam_file%.sam}.bam" 
done

align_with_hisat2()
{
	local R1_FASTQ="$1"
	local R2_FASTQ="$2"
	hisat2 -p "$NUMTHREADS" --mm -1 "$R1_FASTQ" -2 "$R2_FASTQ" -x "$MOUSEREF" --dta
}

SAMPLES=(
	DMSO1_R1_001.fastq.gz
DMSO1_R2_001.fastq.gz
DMSO2_R1_001.fastq.gz
DMSO2_R2_001.fastq.gz
DMSO3_R1_001.fastq.gz
DMSO3_R2_001.fastq.gz
Fingolimod1_R1_001.fastq.gz
Fingolimod1_R2_001.fastq.gz
Fingolimod2_R1_001.fastq.gz
Fingolimod2_R2_001.fastq.gz
Fingolimod3_R1_001.fastq.gz
Fingolimod3_R2_001.fastq.gz
Ozanimod1_R1_001.fastq.gz
Ozanimod1_R2_001.fastq.gz
Ozanimod2_R1_001.fastq.gz
Ozanimod2_R2_001.fastq.gz
Ozanimod3_R1_001.fastq.gz
Ozanimod3_R2_001.fastq.gz
Ponesimod1_R1_001.fastq.gz
Ponesimod1_R2_001.fastq.gz
Ponesimod2_R1_001.fastq.gz
Ponesimod2_R2_001.fastq.gz
Ponesimod3_R1_001.fastq.gz
Ponesimod3_R2_001.fastq.gz
)

cd ~/fastq

for ((i=0; i<${#SAMPLES[@]}; i+=2)); do
	R1_FASTQ="${SAMPLES[$i]}"
	R2_FASTQ="${SAMPLES[$i+1]}"
	align_with_hisat2 "./${R1_FASTQ}" "./${R2_FASTQ}" > "${R1_FASTQ%_R1*}.sam" &
done
