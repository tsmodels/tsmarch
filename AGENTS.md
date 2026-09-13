# Repository conventions

## NEWS.md and issue attribution

When fixing a bug reported through a GitHub issue, always add a `NEWS.md` entry
that cites both the issue number and the GitHub handle of the reporter, so that
the reporter receives attribution. For example:

```
* Fixed <description of the fix> (#7, reported by @kbuchardt).
```

Entries are added under the section for the version currently in `DESCRIPTION`,
creating that section if it does not yet exist. Do not bump the package version
as part of a bug fix unless asked.

## Commit messages

Write commit messages in the repository's existing plain style, describing the
technical change and its rationale. Never mention AI assistance, coding agents
or co-authorship trailers in commits, `NEWS.md`, documentation, vignettes or
code comments.

## Verification

Run these from the package root before considering a change complete:

```
Rscript -e 'devtools::test(".")'
Rscript -e 'devtools::document(".")'
R CMD build . && R CMD check --no-manual tsmarch_<version>.tar.gz
```

`tests/testthat` is the fast suite run by `devtools::test()`. `tests/longtests`
is slow and is not run by default.

Regenerated `man/*.Rd` files are tracked and should be committed with the
roxygen changes that produced them.
