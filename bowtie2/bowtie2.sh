# git@github.com:rdnajac/bowtie2.git
$DATA_DIR=/home/ubuntu/run2/
files=( # {{{
CBP-2-H3K27ac2_S15_L001_R1_001.fastq.gz    
CBP-2-H3K27ac2_S15_L001_R2_001.fastq.gz    
CBP-2-H3K27ac2_S15_L002_R1_001.fastq.gz    
CBP-2-H3K27ac2_S15_L002_R2_001.fastq.gz    
CBP-2-H3K27ac2_S15_L003_R1_001.fastq.gz    
CBP-2-H3K27ac2_S15_L003_R2_001.fastq.gz    
CBP-2-H3K27ac2_S15_L004_R1_001.fastq.gz    
CBP-2-H3K27ac2_S15_L004_R2_001.fastq.gz    
CBP-2-input2_S11_L001_R1_001.fastq.gz      
CBP-2-input2_S11_L001_R2_001.fastq.gz      
CBP-2-input2_S11_L002_R1_001.fastq.gz      
CBP-2-input2_S11_L002_R2_001.fastq.gz      
CBP-2-input2_S11_L003_R1_001.fastq.gz      
CBP-2-input2_S11_L003_R2_001.fastq.gz      
CBP-2-input2_S11_L004_R1_001.fastq.gz      
CBP-2-input2_S11_L004_R2_001.fastq.gz      
Combo-2-H3K27ac2_S16_L001_R1_001.fastq.gz  
Combo-2-H3K27ac2_S16_L001_R2_001.fastq.gz  
Combo-2-H3K27ac2_S16_L002_R1_001.fastq.gz  
Combo-2-H3K27ac2_S16_L002_R2_001.fastq.gz  
Combo-2-H3K27ac2_S16_L003_R1_001.fastq.gz  
Combo-2-H3K27ac2_S16_L003_R2_001.fastq.gz  
Combo-2-H3K27ac2_S16_L004_R1_001.fastq.gz 
Combo-2-H3K27ac2_S16_L004_R2_001.fastq.gz  
Combo-2-input2_S12_L001_R1_001.fastq.gz    
Combo-2-input2_S12_L001_R2_001.fastq.gz    
Combo-2-input2_S12_L002_R1_001.fastq.gz    
Combo-2-input2_S12_L002_R2_001.fastq.gz    
Combo-2-input2_S12_L003_R1_001.fastq.gz    
Combo-2-input2_S12_L003_R2_001.fastq.gz    
Combo-2-input2_S12_L004_R1_001.fastq.gz    
Combo-2-input2_S12_L004_R2_001.fastq.gz    
DMSO-2-H3K27ac2_S13_L001_R1_001.fastq.gz
DMSO-2-H3K27ac2_S13_L001_R2_001.fastq.gz
DMSO-2-H3K27ac2_S13_L002_R1_001.fastq.gz
DMSO-2-H3K27ac2_S13_L002_R2_001.fastq.gz
DMSO-2-H3K27ac2_S13_L003_R1_001.fastq.gz
DMSO-2-H3K27ac2_S13_L003_R2_001.fastq.gz
DMSO-2-H3K27ac2_S13_L004_R1_001.fastq.gz
DMSO-2-H3K27ac2_S13_L004_R2_001.fastq.gz
DMSO-2-input2_S9_L001_R1_001.fastq.gz
DMSO-2-input2_S9_L001_R2_001.fastq.gz
DMSO-2-input2_S9_L002_R1_001.fastq.gz
DMSO-2-input2_S9_L002_R2_001.fastq.gz
DMSO-2-input2_S9_L003_R1_001.fastq.gz
DMSO-2-input2_S9_L003_R2_001.fastq.gz
DMSO-2-input2_S9_L004_R1_001.fastq.gz
DMSO-2-input2_S9_L004_R2_001.fastq.gz
RUX-2-H3K27ac2_S14_L001_R1_001.fastq.gz
RUX-2-H3K27ac2_S14_L001_R2_001.fastq.gz
RUX-2-H3K27ac2_S14_L002_R1_001.fastq.gz
RUX-2-H3K27ac2_S14_L002_R2_001.fastq.gz
RUX-2-H3K27ac2_S14_L003_R1_001.fastq.gz
RUX-2-H3K27ac2_S14_L003_R2_001.fastq.gz
RUX-2-H3K27ac2_S14_L004_R1_001.fastq.gz
RUX-2-H3K27ac2_S14_L004_R2_001.fastq.gz
RUX-2-input2_S10_L001_R1_001.fastq.gz
RUX-2-input2_S10_L001_R1_001_fastqc.zip
RUX-2-input2_S10_L001_R2_001.fastq.gz
RUX-2-input2_S10_L002_R1_001.fastq.gz
RUX-2-input2_S10_L002_R2_001.fastq.gz
RUX-2-input2_S10_L003_R1_001.fastq.gz
RUX-2-input2_S10_L003_R2_001.fastq.gz
RUX-2-input2_S10_L004_R1_001.fastq.gz
RUX-2-input2_S10_L004_R2_001.fastq.gz
) # }}}
$BT2_IDXES=/home/ubuntu/genomes/GCF_000001405.40/bowtie2/hg38


bowtie2 -x $BT2_IDXES -1 $DATA_DIR/$files[0] -2 $DATA_DIR/$files[1]  -S eg.sam
#time bowtie2 -x /home/ubuntu/genomes/GCF_000001405.40/bowtie2/hg38 -1 /home/ubuntu/run2/CBP-2-H3K27ac2_S15_L001_R1_001.fastq.gz -2 /home/ubuntu/run2/CBP-2-H3K27ac2_S15_L001_R2_001.fastq.gz -S eg.sam
#time bowtie2 -x /home/ubuntu/genomes/GCF_000001405.40/bowtie2/hg38 -1 /home/ubuntu/run2/CBP-2-H3K27ac2_S15_L001_R1_001.fastq.gz -2 /home/ubuntu/run2/CBP-2-H3K27ac2_S15_L001_R2_001.fastq.gz -S eg.sam 2> eg.log -p 3


# build a ref genome
bowtie2-build <fasta_file> <basename> -p4
bowtie2-build 


