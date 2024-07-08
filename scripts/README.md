# Scripts

The scripts in this directory are all written in `bash`.

## About

`$ bash --version` prints the version of bash:

```plaintext
GNU bash, version 5.1.16(1)-release (x86_64-pc-linux-gnu)
Copyright (C) 2020 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
```

> [!INFO]
> Being explicit about the shell you are using is a good way
> to ensure future compatibility and maintainability.

### POSIX compatibility

It is safe to assume that we will have at least bash 4.0 available.
Pretend other shells do not exist and write scripts for `bash` only.
These scripts have been validated only on the following machine
(output from `$ uname -a`):

```plaintext
Linux ip-172-31-32-180 6.5.0-1020-aws #20~22.04.1-Ubuntu SMP Wed May  1 16:10:50 UTC 2024 x86_64 x86_64 x86_64 GNU/Linux
```

Writing scripts in `bash` (and being explicit about it) is a good way to
ensure readability and maintainability.

| command | description                    |
| ------- | ------------------------------ |
| `!!`    | repeat last command            |
| `!$`    | last argument of last command  |
| `!^`    | first argument of last command |
| `$_`    | last argument of last command  |
| `$?`    | exit status of last command    |

## Resources

- `man bash`
- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Bash Guide](https://mywiki.wooledge.org/BashGuide)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Safety First!](https://github.com/anordal/shellharden/blob/master/how_to_do_things_safely_in_bash.md)
- [Shell script templates](https://stackoverflow.com/questions/430078/shell-script-templates)

## Extra

- Avoid ["Bashisms"](https://mywiki.wooledge.org/Bashism)
- [arrays](https://mywiki.wooledge.org/BashFAQ/005) are tricky

## Files in this directory

Scripts prefixed with an underscore (`_`) are meant to be thin
wrappers around the binaries they call.

Helper functions ensure the binaries are installed and called
with the correct arguments.

They are meant to explicitly document the exact command
invocation including flags and, when necessary, positional
arguments. This provides a single point of change in case
the default command invocation changes.

They should be as simple as possible and leave error handling
to higher level scripts.

## Tips & Tricks

### Multiline comments

Use `: '` to open and `'` to close.

> [source](https://stackoverflow.com/a/43158193)

```sh
: '
This is a
very neat comment
in bash
'
```

alias for making md5 checksums for all files of a filetype in a directory:

```sh
alias md5sums="find . -type f -oname '*.fastq' -exec md5sum {} + > md5sums.txt"
use ./ instead

function md5sums() {
  md5sum
  find . -type f -name '*.fastq' -exec md5sum {} + > md5sums.txt
}
```
