% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.R
\name{print.tstest.escc}
\alias{print.tstest.escc}
\title{Print method for escc test}
\usage{
\method{print}{tstest.escc}(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  include.decision = FALSE,
  ...
)
}
\arguments{
\item{x}{an object inheriting class \dQuote{tstest.escc}.}

\item{digits}{integer, used for number formatting. Optionally, to avoid
scientific notation, set \sQuote{options(scipen=999)}.}

\item{signif.stars}{logical. If TRUE, \sQuote{significance stars} are printed.}

\item{include.decision}{prints out whether to reject the NULL at the 5\% level
of significance.}

\item{...}{not currently used.}
}
\value{
Invisibly returns the original object.
}
\description{
Print method for objects inheriting class \dQuote{tstest.escc}
}
