# Alignment Tools for Genomic and Transcriptomic Data

| Tool                                                                 | Description                                           | Key Features                                                                                                         | Best For                                        | Source                                                     |
| -------------------------------------------------------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- | ---------------------------------------------------------- |
| [BWA](https://bio-bwa.sourceforge.net/)[^1]                          | Maps short DNA sequences to reference genome          | - Uses Burrows-Wheeler Transform (BWT) for indexing<br>- Efficient for short reads<br>- Supports paired-end reads    | Whole Genome Sequencing (WGS), Exome Sequencing | [GitHub](https://github.com/lh3/bwa)                       |
| [STAR](https://github.com/alexdobin/STAR)[^2]                        | Specialized for RNA-Seq alignment                     | - Uses seed-extension search<br>- Detects novel splice junctions<br>- Fast and accurate for long reads               | RNA-Seq, especially with long reads             | [GitHub](https://github.com/alexdobin/STAR)                |
| [HISAT2](https://daehwankimlab.github.io/hisat2/)[^3]                | Splice-aware aligner for DNA and RNA sequences        | - Uses graph-based alignment<br>- Memory-efficient<br>- Supports both DNA and RNA alignment                          | RNA-Seq, WGS, particularly for large genomes    | [GitHub](https://github.com/DaehwanKimLab/hisat2)          |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)[^4] | Efficient short read aligner                          | - Uses FM-index (similar to BWT)<br>- Supports gapped, local, and paired-end alignment<br>- Memory-efficient         | ChIP-seq, WGS                                   | [GitHub](https://github.com/BenLangmead/bowtie2)           |
| [Subread](https://subread.sourceforge.net/)[^5]                      | Seed-and-vote algorithm-based aligner                 | - Fast and accurate<br>- Supports indel detection<br>- Includes read counting functionality (featureCounts)          | RNA-Seq, DNA-Seq                                | [GitHub](https://github.com/ShiLab-Bioinformatics/subread) |
| [Subjunc](https://subread.sourceforge.net/subjunc.html)[^6]          | Exon-exon junction detector (part of Subread package) | - Detects novel exon-exon junctions<br>- Uses seed-and-vote algorithm<br>- Can be used independently or with Subread | RNA-Seq, specifically for junction detection    | [GitHub](https://github.com/ShiLab-Bioinformatics/subread) |

[^1]: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168](https://pubmed.ncbi.nlm.nih.gov/19451168/)

[^2]: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21. [PMID: 23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)

[^3]: Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12(4):357-360. [PMID: 25751142](https://pubmed.ncbi.nlm.nih.gov/25751142/)

[^4]: Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357-359. [PMID: 22388286](https://pubmed.ncbi.nlm.nih.gov/22388286/)

[^5]: Liao Y, Smyth GK, Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research. 2013;41(10):e108. [PMID: 23558742](https://pubmed.ncbi.nlm.nih.gov/23558742/)

[^6]: Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research. 2019;47(8):e47. [PMID: 30783653](https://pubmed.ncbi.nlm.nih.gov/30783653/)
