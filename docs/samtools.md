# `samtools`

[samtools](http://www.htslib.org/doc/samtools.html) is a suite of programs
for interacting with high-throughput sequencing data that include:

- [Samtools](https:/github.com/samtools/samtools):
  Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
- [BCFtools](https:/github.com/samtools/bcftools):
  Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP
  and short indel sequence variants
- [HTSlib](https:/github.com/samtools/):
  A C library for reading/writing high-throughput sequencing data

## Quick links to the documentation

- [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html):
  sort alignments by leftmost coordinates
- [`samtools view`](https://www.htslib.org/doc/samtools-view.html):
  converts between different formats
- [`samtools flagstat`](https://www.htslib.org/doc/samtools-flagstat.html):
  quickly calculate simple statistics from a BAM file
- [`samtools index`](https://www.htslib.org/doc/samtools-index.html):
  index a BAM file
- [`samtools merge`](https://www.htslib.org/doc/samtools-merge.html):
  merge multiple sorted BAM files
- [`samtools mpileup`](https://www.htslib.org/doc/samtools-mpileup.html):
  multi-way pileup

- TODO: `REF_PATH` and `REF_CACHE`

## `.bam` to `.cram`

CRAM is a compressed version of the BAM format.

- [CRAM format specification](https://samtools.github.io/hts-specs/CRAMv3.pdf)
- [Using samtools to convert BAM to CRAM](https://www.htslib.org/workflow/cram.html)

**_IMPORTANT_**:

1. Alignments should be kept in chromosome/position sort order.
2. The reference must be available at all times.
   Losing it may be equivalent to losing all your read sequences.

> [!TIP]
> when downloading the indexes for pipelines, see the compatibility issue below

```sh

$ samtools faidx download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz
[E::fai_build_core] File truncated at line 1
[E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
[faidx] Could not build fai index download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz.fai

```

decompress with gzip and recompress with bgzip

```sh
gunzip GCA_000001635.9_GRCm39_full_analysis_set.fna.gz && bgzip download/GCA_000001635.9_GRCm39_full_analysis_set.fna
```
