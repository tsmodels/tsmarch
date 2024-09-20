#' Fast combination of n elements, taken m at a time
#'
#' @param x integer vector source for combinations
#' @param m number of elements to choose.
#'
#' @return a matrix with m rows.
#' @details
#' This is a significantly faster version of \code{\link[utils]{combn}}. For
#' an integer vector x of length 1000 and m of 3 which results in a matrix of size
#' 3 x 166,167,000, the speed improvement is about 9 times (11 secs vs 99 secs).
#' This code is written using Armadillo C++ and 64 bit unsingned integers.
#'
#' @export
#'
#' @examples
#' all.equal(combn(1:10, 3), combn_fast(1:10,3))
combn_fast <- function(x, m)
{
    .combn(x, m)
}
