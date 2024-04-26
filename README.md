# cbmf
- my bioinformatic workflows

## software
#!/bin/bash

sudo apt update && sudo apt install -y build-essential
sudo apt-get update 
sudo apt install unzip
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip && sudo ./aws/install && rm -rf awscliv2.zip aws
# add  ~/.aws/config and ~/.aws/credentials

# for fastqc
sudo apt install default-jre
# add bs 


# #!/bin/bash

# Install dependencies for SAMtools and HTSlib
sudo apt update
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev

mkdir src && cd src
# https://github.com/BenLangmead/bowtie2
# https://github.com/rnakato/DROMPAplus.git
# https://github.com/samtools/htslib.git
# https://github.com/samtools/samtools.git 
# https://github.com/samtools/bcftools.git
# https://github.com/madler/zlib.git
# https://github.com/s-andrews/FastQC.git

# Clone and install SAMtools
git clone https://github.com/samtools/samtools.git
cd samtools
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make && sudo make install


# Clone and install HTSlib
git clone https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
autoreconf -i
./configure
make && sudo make install


