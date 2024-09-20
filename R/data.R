#' Dow Jones 30 Constituents Closing Value log Weekly Return
#'
#' Dow Jones 30 Constituents closing value weekly (Friday) log returns
#' from 1987-03-16 to 2009-02-03 from Yahoo Finance. Note that AIG was
#' replaced by KFT (Kraft Foods) on September 22, 2008. This is not reflected
#' in this data set as that would bring the starting date of the data to 2001.
#' When data was not available for a Friday, the closest previous close for
#' which data was available was used.
#'
#' @format ## `dji30retw`
#' A data frame with 1131 rows and 30 columns, with row names being the
#' YYYY-mm-dd date.
#' @source Yahoo Finance
"dji30retw"


#' Global Financial Indices Closing Value log Weekly Return
#'
#' A selection of the log weekly returns (Friday) of 34 global financial indices
#' from Yahoo Finance starting on 1991-12-13 and ending on 2024-06-21.
#'
#' @format ## `globalindices`
#' A data frame with 1698 rows and 34 columns, with row names being the
#' YYYY-mm-dd date.
#' @source Yahoo Finance
"globalindices"
