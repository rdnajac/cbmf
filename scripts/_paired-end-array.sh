#!/bin/bash
#
## A quick and dirty bash library for data struct
## that handles fastq paired-end files using assosiative arrays see in bash 4.0
#!/bin/bash

declare -A fastq_files

populate_fastq_files() {
	local input_folder=$1
	local file base_name pair_name

	for file in "$input_folder"/*.fastq.gz; do
		[[ -f "$file" ]] || continue
		if [[ "$file" =~ _R[12]_ ]]; then
			base_name=$(basename "$file" .fastq.gz)
			pair_name="${base_name%_R[12]}"
			if [[ -n "${fastq_files[$pair_name]}" ]]; then
				if [[ "${fastq_files[$pair_name]}" =~ "$file" ]]; then
					echo "Duplicate found: $file"
				else
					fastq_files[$pair_name]+=" $file"
				fi
			else
				fastq_files[$pair_name]=$file
			fi
		fi
	done
}

check_missing_pairs() {
	local key files
	for key in "${!fastq_files[@]}"; do
		files=("${fastq_files[$key]}")
		if [[ ${#files[@]} -ne 2 ]]; then
			echo "Missing pair for: $key"
		fi
	done
}

print_fastq_files() {
	local key
	for key in "${!fastq_files[@]}"; do
		echo "$key: ${fastq_files[$key]}"
	done
}

# Example usage
# populate_fastq_files "/path/to/your/input_folder"
# check_missing_pairs
# print_fastq_files
