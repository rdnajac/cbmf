#!/bin/bash
function merge_next_improved()
{
	local input_dir="$1"
	cd "$input_dir" || exit
	mkdir -p merged

	for lane in L00{1..4}; do
		for read in R{1..2}; do
			for f in ./*"${lane}"_"$read"_001.fastq.gz; do
				local new_filename="merged/$(basename "$f" | sed "s/L00[1-4]_R[1-2]_001.fastq.gz/${read}.fastq.gz/")"

				echo "Merging $f into $new_filename"
				cat "$f" >>"$new_filename"
			done
		done
	done
}

# for i in 1 2 3; do
# 	(merge_next_improved kalay/run"$i") &
# done

# for i in 1 2 3; do
# 	(cd kalay/run"$i" && md5sum merged/* >merged/md5sum.txt) &
# done
# wait
# cat kalay/run{1..3}/merged/md5sum.txt

for i in 1 2 3; do
	aws s3 cp --recursive kalay/run"$i"/merged s3://lab-aaf-scratch/kalay_run"$i"_merged_fastqs/
done

now we want to set up the folders
~/kalay/merged_fastqs/run{1..5}
mkdir -vp ~/kalay/merged_fastqs/run{1..5}
