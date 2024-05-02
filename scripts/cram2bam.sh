#! /bin/bash
set -euxo pipefail
function convert_format() {
    local input_file=$1
    local output_file=$2
    local fai=$3
    local threads=$(nproc)
    local cmd="samtools view -@ $threads --verbosity=8"
    case "${input_file##*.}" in
        bam)
            [[ "${output_file##*.}" != "cram" ]] && echo "Output must be .cram for BAM inputs." && return 1
            echo "Converting BAM to CRAM..."
            $cmd -C -T $fai $input_file > $output_file
            ;;
        cram)
            [[ "${output_file##*.}" != "bam" ]] && echo "Output must be .bam for CRAM inputs." && return 1
            echo "Converting CRAM to BAM..."
            $cmd -b -t $fai $input_file > $output_file
            ;;
        *)
            echo "Unsupported file type. Please provide a .bam or .cram file."
            return 1
            ;;
    esac

    if [ $? -eq 0 ]; then
        echo "Conversion completed successfully."
    else
        echo "Failed to convert file."
        return 1
    fi
}

for file in DMSO{1..3}.bam; do
  convert_format $file ${file%.bam}.cram ./download/GCA_000001635.9_GRCm39_full_analysis_set.fna
  # write out the full command:
  samtools view -@ $(nproc) --verbosity=8 -C -T download/GCA_000001635.9_GRCm39_full_analysis_set.fna DMSO1.bam > DMS01.cram
done

