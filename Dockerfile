FROM ubuntu:22.04

ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/bin:$PATH

# Install required system packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y \
    bzip2 \
    ca-certificates \
    curl \
    libssl3 \
    libssl-dev \
    libclang-dev \
    libclang-14-dev \
    libxkbcommon-x11-0 \
    && rm -rf /var/lib/apt/lists/*

# Set up R and CRAN repository,
RUN curl -L "https://raw.githubusercontent.com/eddelbuettel/r2u/master/inst/scripts/add_cranapt_jammy.sh" | bash

# Install RStudio Server
RUN curl -LO https://download1.rstudio.org/electron/jammy/amd64/rstudio-2024.04.2-764-amd64.deb \
    && dpkg -i rstudio-2024.04.2-764-amd64.deb && rm rstudio-2024.04.2-764-amd64.deb

# Install micromamba and create conda environment
RUN mkdir -p /opt/micromamba \
    && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /opt/micromamba \
    && /opt/micromamba/bin/micromamba shell init -s bash -p $MAMBA_ROOT_PREFIX \
    && /opt/micromamba/bin/micromamba create -y -n cbmf \
        -c bioconda -c conda-forge \
        bedtools bowtie2 bwa cellbender fastqc gatk4 hisat2 \
        jupyter jupyterlab matplotlib multiqc numpy pandas \
        samtools scipy star stringtie subread \
    && /opt/micromamba/bin/micromamba clean --all --yes

# Expose ports for JupyterLab and RStudio Server
EXPOSE 8888 8787

# Start RStudio Server and JupyterLab
CMD ["/bin/bash", "-c", "source /opt/micromamba/etc/profile.d/mamba.sh 
&& micromamba activate cbmf && service rstudio-server start &&
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root
