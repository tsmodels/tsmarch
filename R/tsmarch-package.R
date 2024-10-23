#' @rawNamespace useDynLib(tsmarch, .registration=TRUE);
#' @keywords internal
#' @import tsmethods
#' @import data.table
#' @import methods
#' @importFrom stats mahalanobis residuals fitted sigma printCoefmat simulate predict AIC BIC approxfun ar cov.wt cov2cor density dnorm ecdf integrate logLik median na.omit optim pnorm qnorm quantile rnorm stepfun vcov
#' @importFrom zoo na.fill coredata index is.zoo `coredata<-` as.zoo
#' @importFrom tsgarch garch_modelspec newsimpact to_multi_estimate
#' @importFrom xts xts as.xts is.xts merge.xts
#' @importFrom sandwich estfun bwNeweyWest vcovHAC vcovOPG bread
#' @importFrom nloptr nloptr
#' @importFrom grDevices hcl.colors heat.colors terrain.colors topo.colors
#' @importFrom Rsolnp solnp
#' @importFrom graphics grid layout lines par axis barplot contour hist image mtext persp points title abline box
#' @importFrom future.apply future_lapply
#' @importFrom future future
#' @importFrom abind abind
#' @importFrom shape drapecol
#' @importFrom numDeriv hessian jacobian grad
#' @importFrom tsdistributions ddist rdist qdist pdist dskewness dkurtosis spd_modelspec pspd dspd qspd rspd nigtransform ghyptransform qgh
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom RcppBessel bessel_k
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom lubridate tz days weeks years `%m+%`
#' @importFrom utils tail head data combn write.table

"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
