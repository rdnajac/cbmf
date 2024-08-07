# ChIP-seq Analysis

## Overview

Chromatin immunoprecipitation sequencing (ChIP-seq) is used to analyze protein
interactions with DNA. It combines chromatin immunoprecipitation (ChIP) with massively
parallel DNA sequencing to identify the binding sites of DNA-associated proteins.

### parse2wig

parse2wig is a tool for converting aligned read files to wiggle format, which is useful for visualization in genome browsers.

### 3. Peak Calling

Peak calling is a crucial step in ChIP-seq analysis, identifying regions of the genome where proteins of interest are bound. Several tools are available:

| Tool                     | Description                                                                                                           | Source                                                       |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| MACS2[^1]                | Model-based Analysis of ChIP-Seq. Captures local biases in sequencing data and models peak shape to improve accuracy. | [GitHub](https://github.com/macs3-project/MACS)              |
| IDR[^2]                  | Irreproducible Discovery Rate. Assesses reproducibility of findings across replicates.                                | [GitHub](https://github.com/nboley/idr)                      |
| phantompeakqualtools[^3] | A package for ChIP-seq quality control metrics, including strand cross-correlation.                                   | [GitHub](https://github.com/kundajelab/phantompeakqualtools) |

### 4. Visualization and Analysis

#### DROMPA

DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a comprehensive tool for easy-to-handle peak calling and visualization in ChIP-seq data analysis[^4].

Key features:

- Peak calling
- Visualization of enrichment profiles
- Annotation of peaks
- Computational analysis and validation of ChIP-seq data

| Feature              | Description                                                           |
| -------------------- | --------------------------------------------------------------------- |
| Peak Calling         | Identifies significant binding sites using a sliding window approach  |
| Visualization        | Generates publication-quality figures of ChIP-seq enrichment profiles |
| Annotation           | Associates peaks with genomic features (e.g., promoters, genes)       |
| Comparative Analysis | Allows comparison of multiple ChIP-seq datasets                       |

## References

[^1]: Zhang Y, Liu T, Meyer CA, et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. [PMID: 18798982](https://pubmed.ncbi.nlm.nih.gov/18798982/)

[^2]: Li Q, Brown JB, Huang H, Bickel PJ. Measuring reproducibility of high-throughput experiments. Ann Appl Stat. 2011;5(3):1752-1779. [DOI: 10.1214/11-AOAS466](https://projecteuclid.org/euclid.aoas/1318514284)

[^3]: Landt SG, Marinov GK, Kundaje A, et al. ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res. 2012;22(9):1813-1831. [PMID: 22955991](https://pubmed.ncbi.nlm.nih.gov/22955991/)

[^4]: Nakato R, Itoh T, Shirahige K. DROMPA: easy-to-handle peak calling and visualization software for the computational analysis and validation of ChIP-seq data. Genes Cells. 2013;18(7):589-601. [PMID: 23672187](https://pubmed.ncbi.nlm.nih.gov/23672187/)
