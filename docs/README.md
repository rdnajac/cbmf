# Doumentation 📝

> Don't document the program; program the document.
>
> — _the internet_

Institute of Electrical and Electronics Engineers (IEEE)
publications on writing effective technical documents:

- [Write a Good Technical Report](https://ieeexplore.ieee.org/document/6448763)
- [Code Documentation](https://ieeexplore.ieee.org/abstract/document/5484109)

## Style

### Writing

The `.md` files in this repository are written using
[GitHub Flavored Markdown \(GFM\)](https://github.github.com/gfm/),
a lightweight markup language with plain text formatting syntax.

> [!NOTE]
> GFM is a superset of the original [Markdown](https://daringfireball.net/projects/markdown/syntax),
> with additional features like tables, task lists, and alerts like this one.

Since the repository is hosted on GitHub, each `README.md` will be rendered as
a webpage when viewing either the folder containing the file or the file itself.[^1]
The idea is to write the files so that they are equally readable as plain text
and when rendered as a webpage on GitHub. For more information, read the guide to
[basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

#### Spell checking with `cspell`

Use [`cspell`](https://cspell.org/) ckage to spell check markdown files:

```sh
npx cspell **/*.md
```

## Code Style

Note to self: stop agonizing over stylistic details and use `prettier`
to automatically and consistently format markdown and html files.

```sh
prettier --write **/*.md **/*.html
```

## Stuff I like

Check out this article on [semantic line breaks](https://sembr.org/).
https://cspell.org/

<!-- Footnotes -->

[^1]: [About READMEs](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes)
[^1]: For more information, read the guide to [basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).
