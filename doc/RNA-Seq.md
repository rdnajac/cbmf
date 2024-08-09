# RNA Sequencing Analysis

RNA sequencing (RNA-seq) is a technique used to analyze the transcriptome -
the complete set of RNA transcripts in a cell. This method involves sequencing
RNA molecules after reverse transcription to cDNA to examine gene expression levels
and identify novel transcripts.

## `featureCounts`

[`featureCounts`](https://subread.sourceforge.net/featureCounts.html)
is a program for counting reads mapped to genomic features, such as genes, exons,and promoters.[^1]
It is part of the [Subread](https://subread.sourceforge.net) package and is widely used for RNA-seq data analysis.

[^1]: Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-30. [PMID: 24227677](https://pubmed.ncbi.nlm.nih.gov/24227677/)

### Input

Aligned reads in BAM format and a GTF file containing genomic features.

Example usage:

```sh
featureCounts -T "$NUM_THREADS" --verbose -t exon -g gene_id --countReadPairs \
		-a "$REF_GTF" -p -P -C -B -o "${OUT_DIR}/${OUT_PREFIX}.tsv" ./*.bam
```

### Output

Feature counts are written to a tab-delimited file (`.tsv` or `.txt`)
with columns for each sample. You can import this file into R or other
statistical software for further analysis.

Feature counts also provide a summary of the number of reads that were
assigned to features, as well as the number of reads that were not assigned
to any feature.

## Tuxedo Suite

The Tuxedo Suite is a collection of tools for transcript-level expression analysis of RNA-seq experiments.[^2]

[^2]: Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. [PMID: 27560171](https://pubmed.ncbi.nlm.nih.gov/27560171/)

| Tool          | Description                                    | Manual                                                                      | Source                                             |
| ------------- | ---------------------------------------------- | --------------------------------------------------------------------------- | -------------------------------------------------- |
| HISAT2[^3]    | Align reads to reference genome                | [manual](https://daehwankimlab.github.io/hisat2/manual/)                    | [source](https://github.com/DaehwanKimLab/hisat2)  |
| StringTie[^4] | Assemble RNA-Seq alignments into transcripts   | [manual](https://ccb.jhu.edu/software/stringtie/index.shtml)                | [source](https://github.com/gpertea/stringtie)     |
| Ballgown[^5]  | Isoform-level differential expression analysis | [manual](https://bioconductor.org/packages/release/bioc/html/ballgown.html) | [source](https://github.com/alyssafrazee/ballgown) |

[^3]: Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12(4):357-360. [PMID: 25751142](https://pubmed.ncbi.nlm.nih.gov/25751142/)

[^4]: Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015;33(3):290-295. [PMID: 25690850](https://pubmed.ncbi.nlm.nih.gov/25690850/)

[^5]: Frazee AC, Pertea G, Jaffe AE, Langmead B, Salzberg SL, Leek JT. Ballgown bridges the gap between transcriptome assembly and expression analysis. Nat Biotechnol. 2015;33(3):243-246. [PMID: 25748911](https://pubmed.ncbi.nlm.nih.gov/25748911/)

> [!TIP]
> StringTie comes packaged with `gffcompare` for comparing and evaluating the
> accuracy of RNA-seq transcript assemblers. Read the
> [manual](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) for details.

### Workflow

1. Align reads to the reference genome and sort to BAM format (HISAT2)
2. Assemble transcripts (StringTie)
3. Prepare for differential expression analysis (Ballgown setup)
4. Perform differential expression analysis
5. Visualize results

Using Ballgown:

1. Filter to remove low-abundance genes
2. Identify differentially expressed transcripts and genes
3. Add gene names
4. Sort by p-value
5. Write results to file

Tips:

- Save the Ballgown object to a file for later use: `save(bg, file="bg.rda")`
- Add log2 fold change to the gene-level data frame:
  `log2fc <- log2(gene$mean2/gene$mean1)`
- Use R packages like ggplot2 or built-in Ballgown functions to create informative visualizations of your results.

## Differential Expression Analysis

Differential expression (DE) involves comparing gene expression levels between different conditions
or groups to identify genes that are up- or down-regulated (i.e., differentially expressed).
