% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ica.R
\name{radical}
\alias{radical}
\title{The Robust Accurate, Direct ICA ALgorithm (RADICAL)}
\usage{
radical(
  X,
  components = NCOL(X),
  demean = TRUE,
  pca_cov = c("ML", "LW", "EWMA"),
  k = 150,
  augment = FALSE,
  replications = 30,
  sigma = 0.175,
  first_eigen = NULL,
  last_eigen = NULL,
  E = NULL,
  D = NULL,
  Z = NULL,
  K = NULL,
  L = NULL,
  seed = NULL,
  trace = FALSE,
  ...
)
}
\arguments{
\item{X}{the data matrix (n x m) where n is the number of samples and m is the number of signals.}

\item{components}{the number of independent components to extract.}

\item{demean}{whether to demean the data.}

\item{pca_cov}{the method used to calculate the covariance matrix for PCA. Options are
\dQuote{ML} (maximum likelihood), \dQuote{LW} (Ledoit-Wolf) shrinkage method and
\dQuote{EWMA} (exponentially weighted moving average).}

\item{k}{the number of angles at which to evaluate the contrast function. The
contrast function will be evaluated at K evenly spaced rotations from -Pi/4 to Pi/4.}

\item{augment}{whether to augment the data (as explained in paper). For large
datasets of >10,000 points this should be set to FALSE.}

\item{replications}{the number of replicated points for each original point if
augment is set to TRUE. The larger the number of points in the data set,
the smaller this value can be. For data sets of 10,000 points or more, point
replication should be de-activated by setting augment to FALSE.}

\item{sigma}{the standard deviation (noise) of the replicated points when using
the augmentation option (which sets the standard deviation for the random normal
generator).}

\item{first_eigen}{This and \code{last_eigen} specify the range for eigenvalues
that are retained, \code{first_eigen} is the index of largest eigenvalue to be
retained. Making use of this option overwrites \code{components}.}

\item{last_eigen}{the index of the last (smallest) eigenvalue to be retained
and overwrites \code{component} argument.}

\item{E}{(optional) Eigen vector from the PCA decomposition. If this
is provided then \code{D} must also be provided.}

\item{D}{(optional) Eigen values from the PCA decomposition.}

\item{Z}{(optional) provided whitened signal matrix. If this is provided then
both \code{K} and \code{L} must also be provided.}

\item{K}{(optional) whitening matrix.}

\item{L}{(optional) de-whitening matrix.}

\item{seed}{the random seed for the random number generator.}

\item{trace}{whether to print out progress information.}

\item{...}{additional arguments to be passed to the covariance matrix calculation.
For arguments passed to the \dQuote{EWMA} method, it optionally takes an additional argument
\code{lambda} which is the exponential decay parameter (default is 0.96).
The \dQuote{LW} takes an additional argument \code{shrink} which is the
shrinkage parameter (default is to calculate this).}
}
\value{
A list with the following components:
\item{A}{the mixing matrix}
\item{W}{the unmixing matrix}
\item{S}{the independent components}
\item{U}{the rotation matrix}
\item{K}{the whitening matrix}
\item{L}{the dewhitening matrix}
\item{C}{the covariance matrix}
\item{Z}{the whitened signal}
\item{mu}{the mean of the mixed signal (X)}
\item{elapsed}{the time taken to run the algorithm}
}
\description{
The ICA algorithm of Learned-Miller (2003), is based on an efficient entropy
estimator (due to Vasicek (1976)) which is robust to outliers and requires
no strong characterization assumptions about the data generating process.
}
\details{
Steps to the general algorithm are as follows (see P.1284 of
Learned-Miller (2003) for specific details of RADICAL implementation):
\enumerate{
\item Demean the data if required: \eqn{M = X - \mu}
\item Calculate the covariance matrix \eqn{\Sigma}  using one of the
methods provided.
\item Use an eigen decomposition to calculate the eigenvalues and eigenvectors
of the covariance matrix: \eqn{\Sigma = E D E'}
\item Select the range of eigenvalues to retain (dimensionality reduction).
\item Calculate the whitening matrix \eqn{K = D^{-1/2}E'} and the dewhitening
matrix \eqn{L = E D^{1/2}}.
\item Whiten the data: \eqn{Z = M K'}. Unwhitening is done by \eqn{M = Z L'}.
\item Run the RADICAL algorithm to calculate the rotation matrix \eqn{U},
the mixing matrix: \eqn{A = U^{-1} L} and the unmixing matrix \eqn{W = K' U}.
\item Calculate the independent components: \eqn{S = M W + \bold{1}\mu W } where
\bold{1} is a matrix of ones with dimension (samples x 1).
}

Notice that in calculating the mixing (A) and unmixing (W) matrices we have combined
the whitening (K) and un-whitening (L) matrices with the rotation matrix \eqn{U}.
}
\note{
Replications carries a significant computational cost. The algorithm
is written in C++ and uses \dQuote{RcppParallel} for the most expensive
calculations. See \code{\link[RcppParallel]{setThreadOptions}} for setting
the number of threads to use.
}
\references{
\insertRef{LearnedMiller2003}{tsmarch}
\insertRef{Ledoit2004}{tsmarch}
\insertRef{Vasicek1976}{tsmarch}
}
