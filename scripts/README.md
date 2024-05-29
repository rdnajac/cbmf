
# Shell scripting

Rules:

1. Always use `#!/bin/bash` as the shebang
2. Prefer naming scripts with `.sh` extension
3. Use `set` so that calling the script with bash does not break it
4. Use `set euxo pipefail` to enable debugging and exit on error
5. Forbid SUID and SGID for security reasons
6. Error messages should be printed to `STDERR` and not `STDOUT`
7. Use `printf` instead of `echo` for portability
8. Start each file with a description of its contents. Every file must have a top-level comment including a brief overview of its contents.

## Bash

- [Bash Reference Manual](https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
- [Bash Best Practices](https://bertvv.github.io/cheat-sheets/Bash.html)
- [Shell script templates](https://stackoverflow.com/questions/430078/shell-script-templates)

