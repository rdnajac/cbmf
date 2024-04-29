function bam2cram() {
    # Usage: bam2cram input.bam output.cram reference.fa

    if [ "$#" -ne 3 ]; then
        echo "Usage: bam2cram <input.bam> <output.cram> <reference.fa>"
        return 1
    fi

    local input_bam=$1
    local output_cram=$2
    local reference=$3

    # Convert BAM to CRAM
    samtools view -C -T $reference -o $output_cram $input_bam

    if [ $? -eq 0 ]; then
        echo "Conversion to CRAM completed successfully."
    else
        echo "Failed to convert BAM to CRAM."
        return 1
    fi
}

