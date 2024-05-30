# Scripts

Bash comes pre-installed...  which version?

## Bash


- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Bash Best Practices](https://bertvv.github.io/cheat-sheets/Bash.html)
- [Shell script templates](https://stackoverflow.com/questions/430078/shell-script-templates)

Guidelines:

1. Always use `#!/bin/bash` as the shebang
2. Prefer naming scripts with `.sh` extension
3. Use `set` so that calling the script with bash does not break it
4. Use `set euxo pipefail` to enable debugging and exit on error
6. Error messages should be printed to `STDERR` and not `STDOUT`
8. Start each file with a description of its contents.

Generally, shell scripts will go here while python and perl scripts go in their respective directories.

## Utils

`multifastqc`

```sh
# run fastqc on all files from an input_dir
fastqc -o "$output_dir" --noextract --memory 1024 -t "$(nproc)" "$input_dir"/*

# zip all files from the output_dir
"$zipped_html_file" "$output_dir"/*.html
```

> [!WARN]
> those variables have to be defined.

SAMtools `conversion_utilities.sh`

``` sh
# aliases for converting sample read files
# `source conversion_utilities.sh` to add them to your environment

alias fastq_to_fasta="sed 'N;x;N;N;x;s/@/>/"
alias paired_to_tab5="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed  -n 'n;h;n;g;N;s/\n/	/g;p' reads_2.fq) > reads_12.tab5"
alias paired_to_tab6="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed 'N;x;N;g;N;s/\n/	/g' reads_2.fq) > reads_12.tab6"
alias paired_to_interleaved="paste -d'\n' <(sed 'N;N;N;s/\n/	/g' reads_1.fq) <(sed 'N;N;N;s/\n/	/g' reads_2.fq) | tr '\t' '\n' > reads_12.fq"
```


Convert between sam and bam formats:
``` sh
samtools view -@ "$MAX_THREADS" -C -T "$FNA" -o "${file%.bam}.cram" "$file"
samtools view -@ "$MAX_THREADS" -b -o "${file%.cram}.bam" "$file"
```

