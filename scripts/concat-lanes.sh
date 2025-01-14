#!/bin/bash
# Function to concatenate FASTQ files and rename them

concatenate_fastq() {
	# Input directory as the first positional parameter
	input_dir="$1"

	# Check if the input directory exists
	if [[ ! -d "$input_dir" ]]; then
		echo "Error: Directory '$input_dir' does not exist."
		exit 1
	fi

	output_dir="./${input_dir}concatenated"

	# Create output directory if it doesn't exist
	mkdir -p "$output_dir"

	# Find all unique prefixes based on the part before _L00X
	find "$input_dir" -name "*.fastq.gz" | sed -E 's/_L00[1-4]_R1_001\.fastq\.gz//' | sort -u | while read -r prefix; do
		echo "Processing prefix: $prefix"

		# Construct output file name by removing `_L00X`
		output_file="$output_dir/$(basename "$prefix")_R1_001.fastq.gz"

		# Concatenate all matching files for this prefix
		find "$input_dir" -name "$(basename "$prefix")_L00[1-4]_R1_001.fastq.gz" -print0 | sort -z | xargs -0 cat > "$output_file"

		# Rename the output file to keep everything up to the first `R#`
		renamed_file=$(echo "$output_file" | sed -E 's/(R[0-9]+).*/\1.fastq.gz/')

		mv "$output_file" "$renamed_file"

		echo "Created and renamed: $renamed_file"
	done
}

# Call the function with the input directory as a parameter
# concatenate_fastq "$1"
concatenate_fastq ./D0-D07
concatenate_fastq ./D0-D14
