# Demultiplexing Illumina sequencing data

## bcl2fastq

> Read the [User Guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf).

If you have raw sequencing data in BCL format, you will need to convert it to
FASTQ format using the bcl2fastq2 Conversion Software.
This step can be skipped if you used Azenta for sequencing,
or if you correctly uploaded a valid sample sheet prior to sequencing.

### Installation

Download bcl2fastq2 Conversion Software v2.20 Installer (Linux rpm) from
[Illumina](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html).

Check out this [post](https://www.biostars.org/p/266897/) for instructions
on how to convert this rpm (Red Hat Package Manager) file
into a deb (Debian Package Manager) file.

```sh
sudo alien -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
```

The `-i` flag installs the package after converting it to a temporary deb file.
