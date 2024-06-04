# Genomes 🧬

A guide to aligning raw reads to a reference genome using the following tools:

- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- STAR

> [!NOTE]
> We only need to use one of these tools to align reads to a reference genome,
> but it is useful to have multiple options available.

## NCBI Datasets and The Genome Reference Consortium

The most recent major releases from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets):

- [human](https://www.ncbi.nlm.nih.gov/grc/human)
- [mouse](https://www.ncbi.nlm.nih.gov/grc/mouse)

Link to the [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)

### Download Data

Download the data from the FTP server:

- [GRCm39 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/)
- [GRCh38 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)

#### mouse

`GCA_000001635.9_GRCm39_full_analysis_set`

| file                     | description            | action     |
| ------------------------ | ---------------------- | ---------- |
| fna.bowtie_index.tar.gz  | Bowtie2 index files    | `tar -xvf` |
| fna.hisat2_index.tar.gz  | HISAT2 index files     | `tar -xvf` |
| fna.fai                  | Samtools index file    | N/A        |
| fna.gz                   | FASTA format sequences | ?          |
| refseq_annotation.gff.gz | GFF3 format annotation | `gunzip`   |
| refseq_annotation.gtf.gz | GTF format annotation  | `gunzip`   |

For more information, refer to [this document](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt)

### Building Indexes

It is not enough to download the reference genome and annotation files.
We also need to build the indexes for the alignment tools.

> [!TIP]
> Use the files contained in the `seqs_for_alignment_pipelines.ucsc_ids` folder
> instead of building the indexes from scratch.

## Stuff you should know

[FAQs - NCBI](https://ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/)

### GCA vs. GCF

> The GenBank (GCA) assembly is an archival record that is owned by the submitter and may or may not include annotation. A RefSeq (GCF) genome assembly represents an NCBI-derived copy of a submitted GenBank (GCA) assembly. RefSeq (GCF) assembly records are maintained by NCBI.[^1]

### GTF vs. GFF

GTF (a specific version of GFF2) and GFF (versions GFF2 and GFF3) are used in
gene annotation, with GFF3 being more advanced. GTF uses key-value pairs for
attributes, while GFF3 allows hierarchical feature relationships.

Both formats have nine fields:

- use `GTF` for simpler annotations
- use `GFF3` for complex annotations (e.g., NCBI, UCSC)
