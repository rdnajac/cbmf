# Python

## Python3.?

## pipx

install pipx to install python packages in isolated environments

apps will not be globally accessible until your PATH is updated.
Run `pipx ensurepath` to automatically add it, or manually modify your PATH in your shell's config file.

## Python Buffering Behavior

> Like other commands, the Python interpreter can produce output that we can then pass to other commands and processes.
>
> However, python uses line buffering only if both stdin and stdout are a console. When we employ file redirection or pipes, Python performs full buffering.[^1]

### Stuff I like

#### ruff

pipx install ruff

[^1]: https://www.baeldung.com/linux/python-buffering-and-tee

```

```
