# Doumentation 📝

[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)

> Don't document the program; program the document.
>
> — _the internet_

The `.md` files in this repository are written using
[GitHub Flavored Markdown \(GFM\)](https://github.github.com/gfm/),
a lightweight markup language with plain text formatting syntax that
is easy to read and write (with the added bonus of automatic rendering on GitHub)[^1].
[^1]: [About READMEs](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes)

> [!NOTE]
> GFM is a superset of the original [Markdown](https://daringfireball.net/projects/markdown/syntax),
> with additional features like tables, footnotes, and alerts like this one.
> For more information, read the guide to [basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

## Code Style

Use `prettier` to automatically and consistently format markdown and html files.

```sh
prettier --write **/*.md **/*.html
```

## Spell checking with `cspell`

Use [`cspell`](https://cspell.org/) to spell check markdown files:

```sh
npx cspell **/*.md
```

## Stuff I like

Check out this article on [semantic line breaks](https://sembr.org/).

https://github.com/hackergrrl/art-of-readme

Institute of Electrical and Electronics Engineers (IEEE)
publications on writing effective technical documents:

- [Write a Good Technical Report](https://ieeexplore.ieee.org/document/6448763)
- [Code Documentation](https://ieeexplore.ieee.org/abstract/document/5484109)
