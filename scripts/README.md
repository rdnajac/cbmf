# Scripts

Forget about POSIX compatibility and any shell thats not `bash`.

It is safe to assume that we will have at least bash 4.0 available.

> [!TIP]
> Writing scripts in `bash` (and being explicit about it) is a good way to
> ensure readability and maintainability.

```sh
$ bash --version
GNU bash, version 5.1.16(1)-release (x86_64-pc-linux-gnu)
```

## Table of Contents

- [About](#about)
- [Resources](#resources)
- [Tips & Tricks](#tips--tricks)

## About

`!#/bin/bash`

## Resources

- `man bash`
- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Bash Guide](https://mywiki.wooledge.org/BashGuide)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Safety First!](https://github.com/anordal/shellharden/blob/master/how_to_do_things_safely_in_bash.md)
- [Shell script templates](https://stackoverflow.com/questions/430078/shell-script-templates)

### Extra

- Avoid ["Bashisms"](https://mywiki.wooledge.org/Bashism)
- [arrays](https://mywiki.wooledge.org/BashFAQ/005) are tricky

## Tips & Tricks

### Multiline comments

Use `: '` to open and `'` to close.

```sh
: '
This is a
very neat comment
in bash
'
```

> [source](https://stackoverflow.com/a/43158193)

## Utils

`multifastqc`

```sh
# run fastqc on all files from an input_dir
fastqc -o "$output_dir" --noextract --memory 1024 -t "$(nproc)" "$input_dir"/*

# zip all files from the output_dir
"$zipped_html_file" "$output_dir"/*.html
```

`conversion_utilities.sh` from SAMtools

```sh
# aliases for converting sample read files
# `source conversion_utilities.sh` to add them to your environment

alias fastq_to_fasta="sed 'N;x;N;N;x;s/@/>/"
alias paired_to_tab5="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed  -n 'n;h;n;g;N;s/\n/	/g;p' reads_2.fq) > reads_12.tab5"
alias paired_to_tab6="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed 'N;x;N;g;N;s/\n/	/g' reads_2.fq) > reads_12.tab6"
alias paired_to_interleaved="paste -d'\n' <(sed 'N;N;N;s/\n/	/g' reads_1.fq) <(sed 'N;N;N;s/\n/	/g' reads_2.fq) | tr '\t' '\n' > reads_12.fq"
```

Convert between sam and bam formats:

```sh
samtools view -@ "$MAX_THREADS" -C -T "$FNA" -o "${file%.bam}.cram" "$file"
samtools view -@ "$MAX_THREADS" -b -o "${file%.cram}.bam" "$file"
# where $FNA is the reference genome in fasta format
```
