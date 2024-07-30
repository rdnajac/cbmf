#!/bin/bash
#
## This script runs featureCounts on the BAM files in the input directory

# mouse
REF_GFF="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gff"
REF_GTF="GCA_000001635.9_GRCm39_full_analysis_set.refseq_annotation.gtf"

# other variavles
NUM_THREADS=$(nproc)

main() {
	cd "$1" || {
		echo "Directory not found."
		exit 1
	}

	featureCounts -T "$NUM_THREADS" --verbose \
	        -t exon -g gene_id --countReadPairs \
		-a "${HOME}/genomes/mouse/${REF_GTF}" \
		-p -P -C -B -o ra_counts.tsv ./*.bam
}

main "$1"

# run R
cat << 'END_COMMENT' | R --no-save
library("DESeq2")
fc <- read.table("ra_counts.tsv", header=TRUE, row.names=1)
colnames(fc) <- gsub("_sorted_markdup.bam", "", colnames(fc))
colnames(fc) <- gsub("^\\.\\.", "", colnames(fc))
# countdata: a table with the fragment counts
# coldata: a table with information about the samples
# trim leading columns
fc <- fc[,6:ncol(fc)]
countdata <- fc
coldata <- colnames(fc)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ cell + dex)
Error in `rownames<-`(`*tmp*`, value = colnames(countData)) : 
  attempt to set 'rownames' on an object with no dimensions



: << 'END_COMMENT'
This is a heredoc (<<) redirected to a NOP command (:).
The single quotes around END_COMMENT are important,
because it disables variable resolving and command resolving
within these lines.  Without the single-quotes around END_COMMENT,
the following two "$()" "$()" commands would get executed:
"$(gibberish command)"
"$(rm -fr mydir)"
comment1
comment2
comment3

source: https://stackoverflow.com/a/46049228

featureCounts: a ultrafast and accurate read summarization program
featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads. It is available in the SourceForge Subread package or the Bioconductor Rsubread package.

Input and output
featureCounts takes as input SAM/BAM files and an annotation file including chromosomal coordinates of features. It outputs numbers of reads assigned to features (or meta-features). It also outputs stat info for the overall summrization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info).

The annotation file should be in either GTF format or a simplified annotation format (SAF) as shown below (columns are tab-delimited):
GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497097	chr1	3660633	3661579	-
...
Features and meta-features
Each entry in the provided annotation file is taken as a feature (e.g. an exon). A meta-feature is the aggregation of a set of features (e.g. a gene). The featureCounts program uses the gene_id attribute available in the GTF format annotation (or the GeneID column in the SAF format annotation) to group features into meta-features, ie. features belonging to the same meta-feature have the same gene identifier.

featureCounts can count reads at either feature level or at meta-feature level. When summarizing reads at meta-feature level, read counts obtained for features included in the same meta-feature will be added up to yield the read count for the corresponding meta-feature.

Overlap between reads and features
A read is said to overlap a feature if at least one read base is found to overlap the feature. For paired-end data, a fragment (or template) is said to overlap a feature if any of the two reads from that fragment is found to overlap the feature.

By default, featureCounts does not count reads overlapping with more than one feature (or more than one meta-feature when summarizing at meta-feature level). Users can use the -O option to instruct featureCounts to count such reads (they will be assigned to all their overlapping features or meta-features).

Note that, when counting at the meta-feature level, reads that overlap multiple features of the same meta-feature are always counted exactly once for that meta-feature, provided there is no overlap with any other meta-feature. For example, an exon-spanning read will be counted only once for the corresponding gene even if it overlaps with more than one exon.

Example commands
Below gives example commands of using featureCounts included in the SourceForge Subread package. For the example commands of using featureCounts in Rsubread package, please see the Subread/Rsubread Users Guide. Not that featureCounts automatically detects the format of input read files (SAM/BAM).

Summarize a single-end read dataset using 5 threads:
featureCounts -T 5 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.sam
Summarize a BAM format dataset:
featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.bam
Summarize multiple datasets at the same time:
featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam
Perform strand-specific read counting (use '-s 2' if reversely stranded):
featureCounts -s 1 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.bam
Summarize paired-end reads and count fragments (instead of reads):
featureCounts -p --countReadPairs -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
Summarize multiple paired-end datasets:
featureCounts -p --countReadPairs -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam
Count the fragments that have fragment length between 50bp and 600bp only:
featureCounts -p --countReadPairs -P -d 50 -D 600 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
Count those fragments that have both ends mapped only:
featureCounts -p --countReadPairs -B -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
Exclude chimeric fragments from fragment counting:
featureCounts -p --countReadPairs -C -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
Citation
Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30.


ubuntu@ip-172-31-32-180:~/cbmf$ featureCounts --help
featureCounts: unrecognized option '--help'

Version 2.0.6

Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ...

## Mandatory arguments:

  -a <string>         Name of an annotation file. GTF/GFF format by default. See
                      -F option for more format information. Inbuilt annotations
                      (SAF format) is available in 'annotation' directory of the
                      package. Gzipped file is also accepted.

  -o <string>         Name of output file including read counts. A separate file
                      including summary statistics of counting results is also
                      included in the output (' \
	\
	end. < string > .summary'). Both files
                      are in tab delimited format.

  input_file1 [input_file2] ...   A list of SAM or BAM format files. They can be
                      either name or location sorted. If no files provided,
                      <stdin> input is expected. Location-sorted paired-end reads
                      are automatically sorted by read names.

## Optional arguments:
# Annotation

  -F <string>         Specify format of the provided annotation file. Acceptable
                      formats include 'GTF' (or compatible GFF format) and
                      'SAF'. 'GTF' by default.  For SAF format, please refer to
                      Users Guide.

  -t <string>         Specify feature type(s) in a GTF annotation. If multiple
                      types are provided, they should be separated by ',' with
                      no space in between. 'exon' by default. Rows in the
                      annotation with a matched feature will be extracted and
                      used for read mapping.

  -g <string>         Specify attribute type in GTF annotation. 'gene_id' by
                      default. Meta-features used for read counting will be
                      extracted from annotation using the provided value.

  --extraAttributes   Extract extra attribute types from the provided GTF
                      annotation and include them in the counting output. These
                      attribute types will not be used to group features. If
                      more than one attribute type is provided they should be
                      separated by comma.

  -A <string>         Provide a chromosome name alias file to match chr names in
                      annotation with those in the reads. This should be a two-
                      column comma-delimited text file. Its first column should
                      include chr names in the annotation and its second column
                      should include chr names in the reads. Chr names are case
                      sensitive. No column header should be included in the
                      file.

# Level of summarization

  -f                  Perform read counting at feature level (eg. counting
                      reads for exons rather than genes).

# Overlap between reads and features

  -O                  Assign reads to all their overlapping meta-features (or
                      features if -f is specified).

  --minOverlap <int>  Minimum number of overlapping bases in a read that is
                      required for read assignment. 1 by default. Number of
                      overlapping bases is counted from both reads if paired
                      end. If a negative value is provided, then a gap of up
                      to specified size will be allowed between read and the
                      feature that the read is assigned to.

  --fracOverlap <float> Minimum fraction of overlapping bases in a read that is
                      required for read assignment. Value should be within range
                      [0,1]. 0 by default. Number of overlapping bases is
                      counted from both reads if paired end. Both this option
                      and '--minOverlap' option need to be satisfied for read
                      assignment.

  --fracOverlapFeature <float> Minimum fraction of overlapping bases in a
                      feature that is required for read assignment. Value
                      should be within range [0,1]. 0 by default.

  --largestOverlap    Assign reads to a meta-feature/feature that has the
                      largest number of overlapping bases.

  --nonOverlap <int>  Maximum number of non-overlapping bases in a read (or a
                      read pair) that is allowed when being assigned to a
                      feature. No limit is set by default.

  --nonOverlapFeature <int> Maximum number of non-overlapping bases in a feature
                      that is allowed in read assignment. No limit is set by
                      default.

  --readExtension5 <int> Reads are extended upstream by <int> bases from their
                      5'

--readExtension3 < int > Reads are extended upstream by from their < int > bases
3' end.

  --read2pos <5:3>    Reduce reads to their 5' most base or 3' most base. Read
                      counting is then performed based on the single base the
                      read is reduced to.

# Multi-mapping reads

  -M                  Multi-mapping reads will also be counted. For a multi-
                      mapping read, all its reported alignments will be
                      counted. The 'NH' tag in BAM/SAM input is used to detect
                      multi-mapping reads.

# Fractional counting

  --fraction          Assign fractional counts to features. This option must
                      be used together with '-M' or '-O' or both. When '-M' is
                      specified, each reported alignment from a multi-mapping
                      read (identified via 'NH' tag) will carry a fractional
                      count of 1/x, instead of 1 (one), where x is the total
                      number of alignments reported for the same read. When '-O'
                      is specified, each overlapping feature will receive a
                      fractional count of 1/y, where y is the total number of
                      features overlapping with the read. When both '-M' and
                      '-O' are specified, each alignment will carry a fractional
                      count of 1/(x*y).

# Read filtering

  -Q <int>            The minimum mapping quality score a read must satisfy in
                      order to be counted. For paired-end reads, at least one
                      end should satisfy this criteria. 0 by default.

  --splitOnly         Count split alignments only (ie. alignments with CIGAR
                      string containing 'N'). An example of split alignments is
                      exon-spanning reads in RNA-seq data.

  --nonSplitOnly      If specified, only non-split alignments (CIGAR strings do
                      not contain letter 'N') will be counted. All the other
                      alignments will be ignored.

  --primary           Count primary alignments only. Primary alignments are
                      identified using bit 0x100 in SAM/BAM FLAG field.

  --ignoreDup         Ignore duplicate reads in read counting. Duplicate reads
                      are identified using bit Ox400 in BAM/SAM FLAG field. The
                      whole read pair is ignored if one of the reads is a
                      duplicate read for paired end data.

# Strandness

  -s <int or string>  Perform strand-specific read counting. A single integer
                      value (applied to all input files) or a string of comma-
                      separated values (applied to each corresponding input
                      file) should be provided. Possible values include:
                      0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                      Default value is 0 (ie. unstranded read counting carried
                      out for all input files).

# Exon-exon junctions

  -J                  Count number of reads supporting each exon-exon junction.
                      Junctions were identified from all the exon-spanning reads
                      in the input (containing 'N' in CIGAR string). Counting
                      results are saved to a file named ' < output_file > .jcounts'

  -G <string>         Provide the name of a FASTA-format file that contains the
                      reference sequences used in read mapping that produced the
                      provided SAM/BAM files. This optional argument can be used
                      with '-J' option to improve read counting for junctions.

# Parameters specific to paired end reads

  -p                  If specified, libraries are assumed to contain paired-end
                      reads. For any library that contains paired-end reads, the
                      'countReadPairs' parameter controls if read pairs or reads
                      should be counted.

  --countReadPairs    If specified, fragments (or templates) will be counted
                      instead of reads. This option is only applicable for
                      paired-end reads. For single-end data, it is ignored.

  -B                  Only count read pairs that have both ends aligned.

  -P                  Check validity of paired-end distance when counting read
                      pairs. Use -d and -D to set thresholds.

  -d <int>            Minimum fragment/template length, 50 by default.

  -D <int>            Maximum fragment/template length, 600 by default.

  -C                  Do not count read pairs that have their two ends mapping
                      to different chromosomes or mapping to same chromosome
                      but on different strands.

  --donotsort         Do not sort reads in BAM/SAM input. Note that reads from
                      the same pair are required to be located next to each
                      other in the input.

# Number of CPU threads

  -T <int>            Number of the threads. 1 by default.

# Read groups

  --byReadGroup       Assign reads by read group. "RG" tag is required to be
                      present in the input BAM/SAM files.

# Long reads

  -L                  Count long reads such as Nanopore and PacBio reads. Long
                      read counting can only run in one thread and only reads
                      (not read-pairs) can be counted. There is no limitation on
                      the number of 'M' operations allowed in a CIGAR string in
                      long read counting.

# Assignment results for each read

  -R <format>         Output detailed assignment results for each read or read-
                      pair. Results are saved to a file that is in one of the
                      following formats: CORE, SAM and BAM. See Users Guide for
                      more info about these formats.

  --Rpath <string>    Specify a directory to save the detailed assignment
                      results. If unspecified, the directory where counting
                      results are saved is used.

# Miscellaneous

  --tmpDir <string>   Directory under which intermediate files are saved (later
                      removed). By default, intermediate files will be saved to
                      the directory specified in '-o' argument.

  --maxMOp <int>      Maximum number of 'M' operations allowed in a CIGAR
                      string. 10 by default. Both 'X' and '=' are treated as 'M'
                      and adjacent 'M' operations are merged in the CIGAR
                      string.

  --verbose           Output verbose information for debugging, such as un-
                      matched chromosome/contig names.

  -v                  Output version of the program.
END_COMMENT
