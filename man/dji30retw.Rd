% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dji30retw}
\alias{dji30retw}
\title{Dow Jones 30 Constituents Closing Value log Weekly Return}
\format{
\subsection{\code{dji30retw}}{

A data frame with 1131 rows and 30 columns, with row names being the
YYYY-mm-dd date.
}
}
\source{
Yahoo Finance
}
\usage{
dji30retw
}
\description{
Dow Jones 30 Constituents closing value weekly (Friday) log returns
from 1987-03-16 to 2009-02-03 from Yahoo Finance. Note that AIG was
replaced by KFT (Kraft Foods) on September 22, 2008. This is not reflected
in this data set as that would bring the starting date of the data to 2001.
When data was not available for a Friday, the closest previous close for
which data was available was used.
}
\keyword{datasets}
