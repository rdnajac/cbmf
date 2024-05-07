# batch processing with samtools

import os
import sys

headers = ["Reads according to samtools", ]

# for all .bam files in the current directory, get the number of reads
def get_read_counts():

    for file in os.listdir(os.getcwd()):
        if file.endswith(".bam"):
            #SAMTOOLS_READS=$(samtools view -c $BAM_FILE)
            value = os.system("samtools view -c " + file)


            os.system("samtools view -c " + file)
        columns[header + " " + str(i)].append(value)


# Extract alignments in SAM format
# echo "Extracting alignments in SAM format..."
# samtools view -h $BAM_FILE > "${BAM_FILE%.bam}.sam"
# Convert BAM to CRAM
# echo "Converting BAM to CRAM..."
# samtools view -C -o "${BAM_FILE%.bam}.cram" $BAM_FILE

# Statistics of the BAM file
echo "Generating statistics..."
samtools flagstat $BAM_FILE


