# Scripts

Generally, shell scripts will go here while python and perl scripts go in their respective directories.

## `conversion_utilities.sh` from samtools

This might be useful:

``` sh
# aliases for converting sample read files
# `source conversion_utilities.sh` to add them to your environment

alias fastq_to_fasta="sed 'N;x;N;N;x;s/@/>/"
alias paired_to_tab5="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed  -n 'n;h;n;g;N;s/\n/	/g;p' reads_2.fq) > reads_12.tab5"
alias paired_to_tab6="paste <(sed 'N;x;N;g;N;s/\n/	/g' reads_1.fq) <(sed 'N;x;N;g;N;s/\n/	/g' reads_2.fq) > reads_12.tab6"
alias paired_to_interleaved="paste -d'\n' <(sed 'N;N;N;s/\n/	/g' reads_1.fq) <(sed 'N;N;N;s/\n/	/g' reads_2.fq) | tr '\t' '\n' > reads_12.fq"
```

## Bash

- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Bash Best Practices](https://bertvv.github.io/cheat-sheets/Bash.html)
- [Shell script templates](https://stackoverflow.com/questions/430078/shell-script-templates)

Rules:

1. Always use `#!/bin/bash` as the shebang
2. Prefer naming scripts with `.sh` extension
3. Use `set` so that calling the script with bash does not break it
4. Use `set euxo pipefail` to enable debugging and exit on error
5. Forbid SUID and SGID for security reasons
6. Error messages should be printed to `STDERR` and not `STDOUT`
7. Use `printf` instead of `echo` for portability
8. Start each file with a description of its contents.
