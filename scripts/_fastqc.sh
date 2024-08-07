#!/bin/bash
#
## Batch processing of FastQC

# Arbitrarily setting the memory limit to 4GB
MEMSIZE=4096

# fastqc attemts to run on every file in the current directory and
# by default, its subdirectories. It will skip any file with an
# extension that is not recognized by fastqc (ie .fastq.gz, .bam)
_fastqc()
{
	local out_dir="${1:-fastqc}"

	# Create output directory
	mkdir -vp "$out_dir"

	## Run FastQC with these flags:
	# -t            Number of threads (pass nproc to use all available cores)
	# --memory 	Memory limit for java processes
	# --noextract	Do not extract the results
	# -o		Output directory (defaults to fastqc)
	fastqc -t "$(nproc)" --memory "$MEMSIZE" --noextract -o "$out_dir" ./*
}

# Run the fastqc function
#
_fastqc "$@"
