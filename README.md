# Table of Contents
1. [Getting Started](#getting-started)
2. [Workflows](#workflows)
    1. [Command Line](#command-line)
3. [Illumina Sequencing](#illumina-sequencing)
    1. [Basespace](#basespace)
    2. [bcl2fastq](#bcl2fastq)
4. [Alignment](#alignment)
    1. [Bowtie2](#bowtie2)



# Getting Started
## prerequisites
``` sh
 sudo apt update && sudo apt upgrade -y
sudo apt install openjdk-11-jdk

```

### Reference genomes from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets) 
|species|assembly|ftp|
|---|---|---|
[human](https://www.ncbi.nlm.nih.gov/grc/human)|GCA_000001405.15_GRCh38 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/
[mouse](https://www.ncbi.nlm.nih.gov/grc/mouse) | GCA_000001635.9_GRCm39 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/

``` sh
# download using ftp
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz

# and build the index files
bowtie2-build GCA_000001405.15_GRCh38_genomic.fna.gz GRCh38
bowtie2-build GCA_000001635.9_GRCm39_genomic.fna.gz GRCm39
```
or

``` sh
# download the prebuilt index files
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index.tar.gz

# and extract the files
tar -xvf GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz
tar -xvf GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index.tar.gz
```

## Amazon Web Services
copy files to/from s3 bucket
``` sh
aws s3 sync <source> <destination>
```
s3 uri format `s3://bucketname/uri`

## bash
- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Bash Best Practices](https://bertvv.github.io/cheat-sheets/Bash.html)
### strict enforcement of error handling
``` sh
  set -o errexit    # abort on nonzero exitstatus
  set -o nounset    # abort on unbound variable
  set -o pipefail   # don't hide errors within pipes

  set -euo pipefail # equivalent
  ```
### debugging 
``` sh 
set -x        # print each command before executing it
```

### environment variables
``` sh
export MOUSEREF=~/genomes/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index
export HUMANREF=~/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index
```
to make them permanent, add them to your `.bashrc`
``` sh
echo 'export MOUSEREF=~/genomes/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index' >> ~/.bashrc
echo 'export HUMANREF=~/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index' >> ~/.bashrc
```

### ssh
`~/.ssh/config`: add the location of the private key file so you don't have to specify it with `-i` each time you connect

### scp
"secure copy" files between computers using `ssh` 
``` sh
# copy file to home dir of remote aws instance
scp localfile.txt aws:~

# copy remote file to pwd
scp aws:~/remotefile.txt .
```

### tmux
Weird stuff can happen with "nested" sessions over `ssh`. If you want to attach to a tmux session on a remote server, you need to use the `-t` flag since `tmux` is not a login shell.
``` sh
ssh aws             # works
ssh aws tmux a      # huh?
ssh aws -t tmux a   # ok  
```

### vim
once you have ssh configured, you can use vim to edit files remotely thanks to the `netrw` plugin that comes shipped with `vim`.  
``` sh
vim scp://aws/remote/path/to/file
```

## Basespace
### Installation
on macOS:
``` sh
# https://github.com/basespace/homebrew-basespace
brew tap basespace/basespace
brew install bs-cli

# then authenticate with your BaseSpace credentials
bs authenticate
```
#### 
|
|---|---|
`accession`|Register a new item
`add`|Add lane or workflow settings
`archive`|Archive entities
`authenticate`|Make an authentication request to BaseSpace (aliases: auth)
`await`|Await completion or status change
`children`|Retrieve child entities (aliases: child)
`clear`|Clear lane or workflow settings
`content`|Show all files in an entity (aliases: contents, dir)
`copy`|Copy an entity to a new location (aliases: duplicate)
`create`|Create an entity
`delete`|Move entities to your trash (aliases: rm, move-to-trash)
`download`|Download files from an entity
`empty`|Empty a store
`export`|Export lane or workflow settings
`get`|Get detailed information about an entity (aliases: info, inspect)
`header`|List headers for an entity (aliases: headers)
`history`|Retrieve account activity records for a user or workgroup
`import`|Import lane or workflow settings
`kill`|Abort entities
`launch`|Execute a task
`link`|Get a direct link for an entity
`list`|List and filter entities (aliases: filter, list-all)
`load`|Load into your environment
`logs`|Retrieve log files (aliases: log)
`rename`|Rename entities
`restore`|Restore items
`revoke`|Invalidate a resource (aliases: expire)
`seqstats`|Sequencing stats
`set`|Set properties in an entity
`translate`|Translate v1 <-> v2 entity IDs
`unarchive`|Restore entities from archival storage
`unlock`|Unlock a locked entity
`update`|Update entities
`upload`|Upload files to an entity
`whoami`|Get information about selected user/configuration
`yield`|Return yield information for an entity (aliases: yields)

### `bs list project`

|Name|Id|
|---|---|
| HEL_H3K27ac_ChIPSeq1 |381715669|
| Untitled from 230306_NS500289_1231_AHW5Y3BGXN |383461079|
| CBP_P300 |388989668| 



``` sh
bs list biosample > biosamples.txt
bs list biosample --project-id 381715669
bs list biosample --project-id 383461079
bs list biosample --project-id 388989668
```

# Alignment
There are a number of tools available for aligning reads to a reference genome. We're using `Bowtie2`.

## Bowtie2
Instead of exporting the path like `export BT2_HOME=/home/ubuntu/src/bowtie2`, the binary is installed in `/usr/local/bin` 

As far as I know, there is no reason not to use the `--mm` flag to use memory-mapped I/O, especially since this reduces overhead with multiple threads. 


### Building Bowtie2
Skip this section if you are not building from source. 

To build bowtie2-build with libsais first make sure that the libsais submodule is available. This can be done in one of the following ways: 
``` sh
# first time cloning
git clone --recursive https://github.com/BenLangmead/bowtie2.git

# existing checkout of bowtie2
git submodule init && git submodule update
```
To build Bowtie2 with SRA and libsais support issue the following command:
``` sh
cmake . -D USE_SRA=1 -D USE_SAIS=1 && cmake --build .
```

### Lambda phage example
from `fastq` (or the compressed `fastq.gz`) to `sam` to `bam` to `sorted.bam`
``` sh
# align example paired-end reads
bowtie2 -x example/index/lambda_virus \
        -1 example/reads/reads_1.fq -2 example/reads/reads_2.fq -S eg2.sam

samtools view -bS eg2.sam > eg2.bam
samtools sort eg2.bam -o eg2.sorted.bam
```
Sorted BAM is a useful format because the alignments are 
- compressed, which is convenient for long-term storage, and 
-  sorted, which is conveneint for variant discovery. 

## Pipeline
instead of using the -S flag with `bowtie2`, just use `>` to redirect the output to a file or pipe it to another program like `samtools` with `|`  

``` sh
bowtie2 -p 8 --mm -x index/lambda_virus \
        -1 reads/reads_1.fq -2 reads/reads_2.fq \
        -p 8 | samtools view -bS - | samtools sort -o sorted.bam
```
or, if we already have the .sam file
``` sh
samtools view -@64 -bS Tet2--RhoaG17V-DMSO-24h-3.sam | samtools sort -@64 -o Tet2--RhoaG17V-DMSO-24h-3.sorted.bam 

# use $nproc threads
samtools view -@ $(nproc) -bS eg2.sam | samtools sort -@ $(nproc) -o eg2.sorted.bam

eg2.sam | samtools sort -o eg2.sorted.bam
```
### optionally...
set an alias for bowtie2
``` sh

# set an alias for bowtie2
alias bt2='bowtie2'

# an alternative to setting the path
alias bowtie2='/usr/local/bin/bowtie2'

# "default" options
alias bowtie2='bowtie2 -p 8 --mm'

# with the reference genome
alias bt2h='bowtie2 -p 8 --mm -x $HUMANREF'
alias bt2m='bowtie2 -p 8 --mm -x $MOUSEREF'

```
##
-p 8 | samtools view -bS - | samtools sort -o sorted.bam


# RNA-Seq
## cufflinks
building from source requires boost libraries
``` sh
# install boost
sudo apt install libboost-all-dev
bjam --toolset=gcc architecture=x86 address_model=64 link=static runtime-link=static stage install
```

#### eigen
a template library for linear algebra
``` sh
# download https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
curl -O https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xvf eigen-3.4.0.tar.gz
ubuntu@ip-172-31-32-180:~$ sudo mv eigen-3.4.0 /usr/local/include/
rm eigen-3.4.0.tar.gz
```

./configure  --with-boost=/usr/include/boost --with-eigen=/usr/local/include/eigen-3.4.0

