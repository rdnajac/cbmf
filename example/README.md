# Examples

## Bowtie2 Lambda phage example

from `fastq` (or the compressed `fastq.gz`) to `sam` to `bam` to `sorted.bam`

```sh
# align example paired-end reads
bowtie2 -x example/index/lambda_virus \
        -1 example/reads/reads_1.fq   \
        -2 example/reads/reads_2.fq | \
        samtools view -bS - |         \
        samtools sort > eg2.sorted.bam
```
