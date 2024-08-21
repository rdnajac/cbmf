align_PE_reads() {
	local aligner="$1" r1_fastq="$2" r2_fastq="$3" reference="$4"
	local align_opts="-p $(nproc) --mm -1 $r1_fastq -2 $r2_fastq -x $reference"
	case "$aligner" in
		hisat2) hisat2 "$align_opts" ;;
		bowtie2) bowtie2 "$align_opts" --dta ;;
		*) err_exit "Invalid aligner: $aligner" ;;
	esac
}

align_to_bam() {
	local aligner="$1" r1_fastq="$2" r2_fastq="$3" reference="$4"
	local bam_file="${r1_fastq%_R1*}_sorted_markdup.bam"
	align_reads "$aligner" "$r1_fastq" "$r2_fastq" "$reference" \
		| samtools sort -n -@ "$(nproc)" - \
		| samtools fixmate -@ "$(nproc)" -m - - \
		| samtools sort -@ "$(nproc)" - \
		| samtools markdup -@ "$(nproc)" - "$bam_file"
}

bam_to_cram() {
	local bam_file="$1"
	local id=$(basename "$bam_file" .bam)
	local cram_file="${id}.cram"
	local reference="/home/ubuntu/genomes/human/GCA_000001405.15_GRCh38_full_analysis_set.fna.bgz"
	samtools view -@"$(nproc)" --cram -T "$reference" "$bam_file" > "$cram_file"
	pr_info "Converted $bam_file to $cram_file"
}

cram_to_bam() {
	local cram_file="$1"
	local id=$(basename "$cram_file" .cram)
	local bam_file="${id}.bam"
	[ -f "$bam_file" ] && mv "$bam_file" "${bam_file}.bak"
	samtools view -@"$(nproc)" --bam "$cram_file" > "$bam_file"
	pr_info "Converted $cram_file to $bam_file"
}

