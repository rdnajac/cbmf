# This script is meant to run STAR using a for loop to iterate through all paired end FASTQs in a given directory
# This script was run using STAR VER 2.7.9a
# Utilized reference genome from GenCode GRCH38p12

# Set variables for scripts
genome="/mnt/data/shares/ref/STAR/GRCh38p12/STAR150"
gtf="/mnt/data/shares/ref/STAR/GRCh38p12/gencode.v31.annotation.gtf"
projectDir="/mnt/data/bobby/scratch/RNAseq/20210809_MYLAHUT78_RomiAfaSyn/30-564761036"
fastqDir="${projectDir}/00_fastq"
alignedDir="${projectDir}/aligned"
nCore=12

cd ${fastqDir}

for i in *_R1_001.fastq.gz; do 
        echo ALIGNING ${i%_R1_001.fastq.gz}
        STAR \
		--runMode alignReads \
        	--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
        	--quantMode GeneCounts \
        	--genomeDir ${genome} \
        	--readFilesIn ${i} ${i%_R1_001.fastq.gz}_R2_001.fastq.gz \
        	--runThreadN 12 \
		--outFileNamePrefix ${alignedDir}/${i%_R1_001.fastq.gz}
        echo ${i%_R1_001.fastq.gz} ALIGNMENT FINISHED
done

STAR --genomeLoad Remove --genomeDir ${genome}
