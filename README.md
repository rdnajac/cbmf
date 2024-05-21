# cbmf 🧬

All the code I don't want to write twice.

> ~~computational biology mad flow~~
>
> combinatorial bioinformatic meta-framework

## Shell scripting

- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Bash Best Practices](https://bertvv.github.io/cheat-sheets/Bash.html)

``` sh
# update and upgrade everything
sudo apt update && sudo apt upgrade -y && sudo apt dist-upgrade -y && sudo apt autoremove -y

# run the last command
!!

# run last arguments of last command
# eg last command was `ls -l /path/to/file`
!$   # /path/to/file
```

### Safety first

Include the following at the top of the script:

``` sh
# strict enforcement of error handling
set -o errexit    # abort on nonzero exit status
set -o nounset    # abort on unbound variable
set -o pipefail   # don't hide errors within pipes

set -euo pipefail # equivalent

# debugging
set -x            # print each command before executing it
```

### Environment Variables:

``` sh
# run in current shell
export MOUSEREF=~/genomes/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index
export HUMANREF=~/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index

# or to make them permanent, add to your `.bashrc`
echo 'export MOUSEREF=~/genomes/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index' >> ~/.bashrc
echo 'export HUMANREF=~/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index' >> ~/.bashrc
```

### Template for script files

Inspiration taken from posts here: [Shell script templates](https://stackoverflow.com/questions/430078/shell-script-templates)

``` sh
#!/bin/bash
## Description:
## Author: Ryan D. Najac
## Last modified: $(date +"%Y-%m-%d")

# Enforce strict error handling and print each command
set -euxo pipefail

# Get the script's filename for usage and debugging
SCRIPT_NAME=$(basename "$BASH_SOURCE" .sh)

# Exit function with error message and status code
bail() {
    echo -ne "$1" >&2
    exit ${2:-1}
}

HELP_MSG="Usage: $SCRIPT_NAME <working_directory>\n
Options:
  -h    Display this help and exit
Example:
  $SCRIPT_NAME /path/to/data   # Run summary and cleaning scripts in /path/to/data
"

usage() {
    local status=2
    if [ "$1" -eq "$1" ] 2>/dev/null; then
        status=$1
        shift
    fi
    bail "${1}${HELP_MSG}" $status
}

# Parse command-line options
while getopts "h" opt; do
    case $opt in
        h) usage 0 ;;
        ?) usage "Invalid option: -$OPTARG \n" ;;
    esac
done

# Validate the number of arguments
shift $((OPTIND - 1)) && [[ $# -lt 1 ]] && usage "too few arguments\n"
```

### Functions

Use '\$' to refer to positional arguments

``` sh
cp1ec2() { scp aws:~/${1} ~/Downloads/ }
cp2ec2() { scp "${1}" aws:~/${2:-.} }
sync-from-ec2() { rsync -avz --progress aws:~/path/to/dir/ . }
```

## Remote access

- Connect to and interact the ec2 instance using `ssh` and `scp`
- Bonus points if you `netrw` to edit files remotely with `vim`
- TODO `tmux` for persistent sessions
- TODO `rsync` for keeping directories in sync
- TODO add section on ftp and downloading from Azenta

---
### `aws`
- [AWS CLI Command Reference](https://docs.aws.amazon.com/cli/latest/index.html)
- [AWS S3 CLI Command Reference](https://docs.aws.amazon.com/cli/latest/reference/s3/index.html)
- s3 uri format `s3://bucketname/uri`
``` sh
# copy files to/from s3 bucket
aws s3 sync <source> <destination>

# copy all files of a certain type to local
aws s3 cp s3://bucketname/uri/*.fastq.gz .
```
---
### `ssh`
edit `~/.ssh/config`
- add the location of the private key file so you don't have to specify it with `-i` each time you connect
- add the `User` and `Hostname` fields so you don't have to specify them each time
``` sh
# example ~/.ssh/config
Host aws
    Hostname ec2-3-235-0-0.compute-1.amazonaws.com
    User ubuntu
    IdentityFile ~/.ssh/aws.pem
```
---
### `scp`
- "secure copy" files between computers using `ssh`
- this means the config from `~/.ssh/config` is used
``` sh
# copy file to home dir of remote aws instance
scp localfile.txt aws:~

# copy remote file to pwd
scp aws:~/remotefile.txt .
```
---

### `rsync`
- useful for keeping all files in a directory up-to-date
``` sh
# copy all files in a directory to a remote server
rsync -avz --progress /path/to/local/dir/ aws:/path/to/remote/dir/
```

### `tmux`
Weird stuff can happen with "nested" sessions over `ssh`. If you want to attach to a tmux session on a remote server, you need to use the `-t` flag since `tmux` is not a login shell.
``` sh
ssh aws             # works
ssh aws tmux a      # huh?
ssh aws -t tmux a   # ok
```

#### `vim`

Once you have ssh configured, you can use vim to edit files remotely thanks to the `netrw` plugin that comes shipped with `vim`.
``` sh
vim scp://aws/remote/path/to/file
Copy current vim buffer to remote server
``` vim
:!scp % aws:~/path/to/remote/file
```
---

## Alignment Workflow

## Reference Genomes

### `samtools`

> Samtools is a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:
> * [Samtools](https:/github.com/samtools/samtools): Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
> * [BCFtools](https:/github.com/samtools/bcftools): Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants
> * [HTSlib](https:/github.com/samtools/): A C library for reading/writing high-throughput sequencing data

#### key commands
with links to the documentation
- [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html) - sort alignments by leftmost coordinates
- [`samtools view`](https://www.htslib.org/doc/samtools-view.html) - converts between different formats
- [`samtools flagstat`](https://www.htslib.org/doc/samtools-flagstat.html) - quickly calculate simple statistics from a BAM file
- [`samtools index`](https://www.htslib.org/doc/samtools-index.html) - index a BAM file
- [`samtools merge`](https://www.htslib.org/doc/samtools-merge.html) - merge multiple sorted BAM files
- [`samtools mpileup`](https://www.htslib.org/doc/samtools-mpileup.html) - multi-way pileup

TODO merge and test QC

### FastQC
FastQC is a quality control tool for high throughput sequence data. It reads in sequence data in a variety of formats and can either provide an interactive application to review the results or create an HTML report.

``` sh
# run fastqc on a fastq file
fastqc -o outputdir/ inputfile.fastq
```
FastQC runs on Java
``` sh
# on ubuntu
sudo apt install openjdk-11-jdk
```

#### Automation
- script to run `fastqc` on all `fastq` files in a directory
- script to parse the output of `fastqc` and generate a summary report

##### `multifastqc.sh`
run `fastqc` on all `fastq` files in a directory

##### `qcsummary.py`

### Reference genomes
Arbitrainly using the most recent major releases of the human and mouse genomes  from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets)
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
optionally, set environment variables
``` sh
export MOUSEREF=~/genomes/GCA_000001635.9_GRCm39_full_analysis_set.fna.bowtie_index
export HUMANREF=~/genomes/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index
```

### Hisat2





Usage: 

### Bowtie2
Instead of exporting the path like `export BT2_HOME=/home/ubuntu/src/bowtie2`, install the binary by `sudo cp`ing it to `/usr/local/bin`

> **Note**: Use the `--mm` flag to use memory-mapped I/O, especially since this reduces overhead with multiple threads.
>
> see: `man mmap`

#### Building Bowtie2
Skip this section if you are not building from source.

Building from source gives the options for:
1. basic automatic dependency management and static linkage of `zstd` and `zlib`
    - `make static-libs && make STATIC_BUILD=1`
2. SRA (Sequence Read Archive) support
    - `make sra-deps && make USE_SRA=1.`
3. libsais support: a state-of-the-art suffix array construction algorithm that speeds up the building of the index
    - `[g]make libsais USE_SAIS_OPENMP=1 *`
    - this is important so it can be multithreaded

To build bowtie2-build with libsais first make sure that the libsais submodule is available. This can be done in one of the following ways:
``` sh
# first time cloning
git clone --recursive https://github.com/BenLangmead/bowtie2.git

# existing checkout of bowtie2
git submodule init && git submodule update
```

#### Building with CMake
To build Bowtie2 with SRA and libsais support:
``` sh
cmake . -D USE_SRA=1 -D USE_SAIS=1 && cmake --build .
```

#### Lambda phage example
from `fastq` (or the compressed `fastq.gz`) to `sam` to `bam` to `sorted.bam`
``` sh
# align example paired-end reads
bowtie2 -x example/index/lambda_virus \
        -1 example/reads/reads_1.fq   \
        -2 example/reads/reads_2.fq | \
        samtools view -bS - |         \
        samtools sort > eg2.sorted.bam
```
Sorted BAM is a useful format because the alignments are
- compressed, which is convenient for long-term storage, and
- sorted, which is conveneint for variant discovery.

---
### `samtools`
- [samtools](http://www.htslib.org/doc/samtools.html) is a suite of programs for interacting with high-throughput sequencing data.
- TODO: `REF_PATH` and `REF_CACHE`

### `.bam` to `.cram`
CRAM is a compressed version of the BAM format that is more efficient for long-term storage. It is a good idea to convert BAM files to CRAM files for long-term storage.
- [CRAM format specification](https://samtools.github.io/hts-specs/CRAMv3.pdf)
- [Using samtools to convert BAM to CRAM](https://www.htslib.org/workflow/cram.html)

***IMPORTANT***:
1. Alignments should be kept in chromosome/position sort order.
2. The reference must be available at all times. Losing it may be equivalent to losing all your read sequences.

### `scripts/`
Work in progress...

## RNA-seq: work in progress...
Thus far, the pipeline has been largely protocol agnostic. The next steps are to incorporate the tools necessary for RNA-seq analysis.
### `cufflinks`
cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples.
Building from source requires:
- [boost libraries](https://www.boost.org/)
- [eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html)
``` sh
sudo apt install libboost-all-dev

# bjam is the build system for boost
bjam --toolset=gcc architecture=x86 address_model=64 link=static runtime-link=static stage install
```
#### eigen
a template library for linear algebra
``` sh
curl -O https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xvf eigen-3.4.0.tar.gz
sudo mv eigen-3.4.0 /usr/local/include  # TODO: check this
rm -f eigen-3.4.0.tar.gz

./configure  --with-boost=/usr/include/boost --with-eigen=/usr/local/include/eigen-3.4.0
```

## Data analysis
### [pandas](https://pandas.pydata.org/)
``` sh
# ubuntu install
sudo apt install python3-pandas
# ideally, use a virtual environment
# TODO: add venv setup
```
### Internal QC metrics
- Sample ID
- Yield (Mbases)
- Mean Quality Score
- % Bases >= 30
- [paired-end] Reads
- Mapped Reads
- % Mapped Reads
- Duplicate Reads
- % Duplicate Reads (of all reads)
- Unique Reads
- % Unique Reads (of all reads)
- % Unique (of mapped)
-

## `import re`
- [regex101](https://regex101.com/)
- Python's `re` module is used for regular expressions
example to match sample IDs from a list of files
``` python
re_fcontent = re.compile('Anc-(?P<sample>\d+)-(?P<rep>[AB])_(?P<GE>S\d+)_(?P<lane>L\d\d\d)_(?P<direction>R\d)_fastqc')
```

## Extra
- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

### md5

### building indexin samtools
when downloading the indexes for pipelines, see the compatibility issue below
```

ubuntu@ip-172-31-32-180:~/rnaseq/bam2cram$ samtools !!
samtools faidx download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz
[E::fai_build_core] File truncated at line 1
[E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
[faidx] Could not build fai index download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz.fai
```

decompress with gzip and recompress with bgzip
```
gunzip GCA_000001635.9_GRCm39_full_analysis_set.fna.gz && bgzip download/GCA_000001635.9_GRCm39_full_analysis_set.fna


