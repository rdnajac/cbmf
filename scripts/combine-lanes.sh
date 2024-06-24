#!/bin/bash
#
# This script will combine all fastq files from the same sample
# that are in different lanes (nextseq runs) into a single file.

# Usage: ./combine-lanes.sh /path/to/samples
cd "$1" || exit 1

# First, let's create a directory to store the merged files
mkdir -p merged

for r11 in ./*L001_R1_001.fastq.gz; do cat "$r11" merged/"${r11/L001_/}"; done
for r12 in ./*L002_R1_001.fastq.gz; do cat "$r12" >>merged/"${r12/L002_/}" & done
for r13 in ./*L003_R1_001.fastq.gz; do cat "$r13" >>merged/"${r13/L003_/}" & done
for r14 in ./*L004_R1_001.fastq.gz; do cat "$r14" >>merged/"${r14/L004_/}" & done

for r21 in ./*L001_R2_001.fastq.gz; do cp "$r21" merged/"${r21/L001_/}"; done
for r22 in ./*L002_R2_001.fastq.gz; do cat "$r22" >>merged/"${r22/L002_/}" & done
for r23 in ./*L003_R2_001.fastq.gz; do cat "$r23" >>merged/"${r23/L003_/}" & done
for r24 in ./*L004_R2_001.fastq.gz; do cat "$r24" >>merged/"${r24/L004_/}" & done

# do it agin but so they run in the bg
for r1 in ./*L002_R1_001.fastq.gz; do cat "$r1" >>merged/"${r1/L002_/}" & done

for r2 in ./*L002_R2_001.fastq.gz; do cat "$r2" >>merged/"${r2/L002_/}" & done
for r1 in ./*L003_R1_001.fastq.gz; do cat "$r1" >>merged/"${r1/L003_/}"; done
for r2 in ./*L003_R2_001.fastq.gz; do cat "$r2" >>merged/"${r2/L003_/}"; done
for r1 in ./*L004_R1_001.fastq.gz; do cat "$r1" >>merged/"${r1/L004_/}"; done
for r2 in ./*L004_R2_001.fastq.gz; do cat "$r2" >>merged/"${r2/L004_/}"; done

for r1 in ./*L00{2..4}_R1_001.fastq.gz; do echo cat "$r1" ./*"${r1/L00[2-4]_/}"; done
# no, cat it to the file in the merged folder
#
for f in ./*L00{2..4}_R1_001.fastq.gz; do echo cat "$f" ./*"${f/L00[2-4]_/}" >merged/"${f/L00[2-4]_/}"; done

# finally...
for r11 in ./*L001_R1_001.fastq.gz; do cat "$r12"  >merged/"${r12/L002_/}" & done && wait
for r12 in ./*L002_R1_001.fastq.gz; do cat "$r12" >>merged/"${r12/L002_/}" & done && wait
for r13 in ./*L003_R1_001.fastq.gz; do cat "$r13" >>merged/"${r13/L003_/}" & done && wait
for r14 in ./*L004_R1_001.fastq.gz; do cat "$r14" >>merged/"${r14/L004_/}" & done &&wait

# Lets make a function to merge 4 lanes

# again but suffixes are alwaus _001.fastq.gz
function merge_nextseq_lanes()
{
	local merged_dir="./merged"
	mkdir -p "$merged_dir" || exit 1

	for lane in "
		local r1=""${lane}_R1_001.fastq.gz
		local r2=""${lane}_R2_001.fastq.gz
		cp ""$r1 ""$merged_dir/${r1/_L00/_}
		cat ""$r2 >>""$merged_dir/${r2/_L00/_} &
	done
	wait
}

function merge_next_improved()
{
	local input_dir=""$1
	cd ""$input_dir || exit
	mkdir -p merged

	for lane in L00{1..4}; do
		for read in R{1..2}; do
			for f in ./*""${lane}_""${read}_001.fastq.gz; do
				local new_filename="merged/"$(basename "$f" | sed "s/L00[1-4]_R[1-2]_001.fastq.gz/${read}.fastq.gz/")

				echo "Merging "$f" into "$new_filename
				cat ""$f >>""$new_filename
			done
		done
	done
}

# for i in 1 2 3; do
# 	(merge_next_improved kalay/run""$i) &
# done

# for i in 1 2 3; do
# 	(cd kalay/run""$i && md5sum merged/* >merged/md5sum.txt) &
# done
# wait
# cat kalay/run{1..3}/merged/md5sum.txt

for i in 1 2 3; do
	aws s3 cp --recursive kalay/run""$i/merged s3://lab-aaf-scratch/kalay_run""${i}_merged_fastqs/
done

now we want to set up the folders
~/kalay/merged_fastqs/run{1..5}
mkdir -vp ~/kalay/merged_fastqs/run{1..5}
