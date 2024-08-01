#!/bin/bash
#!/bin/bash
#
## A quick and dirty bash library for handling fastq paired-end files using associative arrays in bash 4.0+
set -euo pipefail
# set -x

declare -A fastq_files

change_dir_or_exit() {
	local dir=$1
	cd "$dir" || exit 1
}

populate_fastq_files() {
	for file in ./*.fastq.gz; do
		if [[ "$file" =~ _R[12]\.fastq\.gz ]]; then
			base_name=$(basename "$file" .fastq.gz)
			key="${base_name%_R[12]}"
			if [[ "$base_name" =~ _R1 ]]; then
				fastq_files["$key"]="${file}"
			elif [[ "$base_name" =~ _R2 ]]; then
				if [[ -n "${fastq_files["$key"]}" ]]; then
					fastq_files["$key"]="${fastq_files["$key"]} ${file}"
				else
					fastq_files["$key"]="${file}"
				fi
			fi
		fi
	done
}

check_missing_pairs() {
	local key
	for key in "${!fastq_files[@]}"; do
		files="${fastq_files[$key]}"
		if [[ ! "$files" =~ " " ]]; then
			echo "Missing pair for: $key"
			exit 2
		fi
	done
}

print_fastq_files() {
	printf '%s\n' "${fastq_files[@]}"
}
# main
