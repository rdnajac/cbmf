# Genomes

A guide to aligning reads to a reference genome using the following tools:

- Hisat2
- Bowtie2
- STAR

Using the most recent major releases of the human and mouse genomes from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets)

- [human](https://www.ncbi.nlm.nih.gov/grc/human)
- [mouse](https://www.ncbi.nlm.nih.gov/grc/mouse)

## Download Data

The Genome Reference Consortium (GRC) provides the latest major releases of the human and mouse genomes at the following FTP sites:

- [GRCm39 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/)
- [GRCh38 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)

Use the files contained in the `seqs_for_alignment_pipelines.ucsc_ids` to skip the building of the index files.

### GCA vs. GCF

> The GenBank (GCA) assembly is an archival record that is owned by the submitter and may or may not include annotation. A RefSeq (GCF) genome assembly represents an NCBI-derived copy of a submitted GenBank (GCA) assembly. RefSeq (GCF) assembly records are maintained by NCBI.[^1]

### `seqs_for_alignment_pipelines.ucsc_ids`

#### mouse

GCA_000001635.9_GRCm39_full_analysis_set.

| file                     | description            | action     |
| ------------------------ | ---------------------- | ---------- |
| fna.bowtie_index.tar.gz  | Bowtie2 index files    | `tar -xvf` |
| fna.hisat2_index.tar.gz  | HISAT2 index files     | `tar -xvf` |
| fna.fai                  | Samtools index file    | N/A        |
| fna.gz                   | FASTA format sequences | ?          |
| refseq_annotation.gff.gz | GFF3 format annotation | `gunzip`   |
| refseq_annotation.gtf.gz | GTF format annotation  | `gunzip`   |

For more information, refer to [this document](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt)

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

```sh
# first time cloning
git clone --recursive https://github.com/BenLangmead/bowtie2.git

# existing checkout of bowtie2
git submodule init && git submodule update
```

#### Building with CMake

To build Bowtie2 with SRA and libsais support:

```sh
cmake . -D USE_SRA=1 -D USE_SAIS=1 && cmake --build .
```

Sorted BAM is a useful format because the alignments are

- compressed, which is convenient for long-term storage, and
- sorted, which is conveneint for variant discovery.

<!-- footnotes -->

[^1]: [FAQs - NCBI](<https://www.ncbinlm.nih.gov/datasets/docs/v2/troubleshooting/faq/#:~:text=The%20GenBank%20(GCA)%20assembly%20is,records%20are%20maintained%20by%20NCBI>)
