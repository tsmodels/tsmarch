#' The Robust Accurate, Direct ICA ALgorithm (RADICAL)
#' @description
#' The ICA algorithm of Learned-Miller (2003), is based on an efficient entropy
#' estimator (due to Vasicek (1976)) which is robust to outliers and requires
#' no strong characterization assumptions about the data generating process.
#' @param X the data matrix (n x m) where n is the number of samples and m is the number of signals.
#' @param components the number of independent components to extract.
#' @param demean whether to demean the data.
#' @param pca_cov the method used to calculate the covariance matrix for PCA. Options are
#' \dQuote{ML} (maximum likelihood), \dQuote{LW} (Ledoit-Wolf) shrinkage method and
#' \dQuote{EWMA} (exponentially weighted moving average).
#' @param k the number of angles at which to evaluate the contrast function. The
#' contrast function will be evaluated at K evenly spaced rotations from -Pi/4 to Pi/4.
#' @param augment whether to augment the data (as explained in paper). For large
#' datasets of >10,000 points this should be set to FALSE.
#' @param replications the number of replicated points for each original point if
#' augment is set to TRUE. The larger the number of points in the data set,
#' the smaller this value can be. For data sets of 10,000 points or more, point
#' replication should be de-activated by setting augment to FALSE.
#' @param sigma the standard deviation (noise) of the replicated points when using
#' the augmentation option (which sets the standard deviation for the random normal
#' generator).
#' @param first_eigen This and \code{last_eigen} specify the range for eigenvalues
#' that are retained, \code{first_eigen} is the index of largest eigenvalue to be
#' retained. Making use of this option overwrites \code{components}.
#' @param last_eigen the index of the last (smallest) eigenvalue to be retained
#' and overwrites \code{component} argument.
#' @param E (optional) Eigen vector from the PCA decomposition. If this
#' is provided then \code{D} must also be provided.
#' @param D (optional) Eigen values from the PCA decomposition.
#' @param Z (optional) provided whitened signal matrix. If this is provided then
#' both \code{K} and \code{L} must also be provided.
#' @param K (optional) whitening matrix.
#' @param L (optional) de-whitening matrix.
#' @param seed the random seed for the random number generator.
#' @param trace whether to print out progress information.
#' @param ... additional arguments to be passed to the covariance matrix calculation.
#' For arguments passed to the \dQuote{EWMA} method, it optionally takes an additional argument
#' \code{lambda} which is the exponential decay parameter (default is 0.96).
#' The \dQuote{LW} takes an additional argument \code{shrink} which is the
#' shrinkage parameter (default is to calculate this).
#' @details
#' Steps to the general algorithm are as follows (see P.1284 of
#' Learned-Miller (2003) for specific details of RADICAL implementation):
#' \enumerate{
#' \item Demean the data if required: \eqn{M = X - \mu}
#' \item Calculate the covariance matrix \eqn{\Sigma}  using one of the
#' methods provided.
#' \item Use an eigen decomposition to calculate the eigenvalues and eigenvectors
#' of the covariance matrix: \eqn{\Sigma = E D E'}
#' \item Select the range of eigenvalues to retain (dimensionality reduction).
#' \item Calculate the whitening matrix \eqn{K = D^{-1/2}E'} and the dewhitening
#' matrix \eqn{L = E D^{1/2}}.
#' \item Whiten the data: \eqn{Z = M K'}. Unwhitening is done by \eqn{M = Z L'}.
#' \item Run the RADICAL algorithm to calculate the rotation matrix \eqn{U},
#' the mixing matrix: \eqn{A = U^{-1} L} and the unmixing matrix \eqn{W = K' U}.
#' \item Calculate the independent components: \eqn{S = M W + \bold{1}\mu W } where
#' \bold{1} is a matrix of ones with dimension (samples x 1).
#' }
#'
#' Notice that in calculating the mixing (A) and unmixing (W) matrices we have combined
#' the whitening (K) and un-whitening (L) matrices with the rotation matrix \eqn{U}.
#' @note
#' Replications carries a significant computational cost. The algorithm
#' is written in C++ and uses \dQuote{RcppParallel} for the most expensive
#' calculations. See \code{\link[RcppParallel]{setThreadOptions}} for setting
#' the number of threads to use.
#' @references
#' \insertRef{LearnedMiller2003}{tsmarch}
#' \insertRef{Ledoit2004}{tsmarch}
#' \insertRef{Vasicek1976}{tsmarch}
#' @returns A list with the following components:
#' \item{A}{the mixing matrix}
#' \item{W}{the unmixing matrix}
#' \item{S}{the independent components}
#' \item{U}{the rotation matrix}
#' \item{K}{the whitening matrix}
#' \item{L}{the dewhitening matrix}
#' \item{C}{the covariance matrix}
#' \item{Z}{the whitened signal}
#' \item{mu}{the mean of the mixed signal (X)}
#' \item{elapsed}{the time taken to run the algorithm}
#' @export
#'
radical <- function(X, components = NCOL(X), demean = TRUE, pca_cov = c("ML", "LW", "EWMA"),
                    k = 150, augment = FALSE, replications = 30, sigma = 0.175,
                    first_eigen = NULL, last_eigen = NULL, E = NULL, D = NULL,
                    Z = NULL, K = NULL, L = NULL, seed = NULL, trace = FALSE, ...)
{
    # A: mixing matrix
    # W: unmixing matrix
    # S: independent components
    # K: whitening matrix
    # L: dewhitening matrix
    # Z: whitened signal
    # M: mixed_signal
    # mu: mixed_mean
    # C: covariance matrix

    tic <- Sys.time()
    if (!is.null(seed)) set.seed(seed)

    m <- NCOL(X)
    n <- NROW(X)
    if (m > n && trace) warning("\nsignal matrix (X) may be orientated in the wrong way.")
    pca_cov <- match.arg(pca_cov[1], c("ML", "LW", "EWMA"))
    if (demean) {
        # don't use colMeans...losses accuracy
        mu <- apply(X, 2, mean)
        M <- sweep(X, 2, mu, FUN = "-")
    } else {
        mu <- rep(0, m)
        M <- X
    }
    if (components > m) stop("\ncannot choose more components than signals.")
    if (is.null(first_eigen)) first_eigen <- 1
    if (is.null(last_eigen)) last_eigen <- components
    if (!is.null(first_eigen)) {
        if (is.null(last_eigen)) {
            last_eigen <- last_eigen + components
            if (last_eigen > m) stop("\nfirst_eigen + components is greater than number of signals.")
        } else {
            if (last_eigen > m) stop("\nlast_eigen greater than number of signals.")
            components <- last_eigen - first_eigen + 1
            if (first_eigen < 1) stop("\ncannot have zero factors in PCA.")
        }
    }

    # PCA E (m x f)
    if (!is.null(E)) {
        if (is.null(D)) stop("\nD required when passing E.")
        if (!is.matrix(E)) stop("\nE must be an m (signals) x f (components) matrix.")
        if (NROW(E) != m) stop("\nE must be an m (signals) x f (components) matrix.")
        if (NCOL(E) != components) stop("\nE must be an m (signals) x f (components) matrix.")
    }

    # PCA D (f by f)
    if (!is.null(D)) {
        if (is.null(D)) stop("\nE required when passing D.")
        if (!is.matrix(D)) stop("\nD must be an f (components) by f (components) diagonal matrix.")
        if (NROW(D) != components) stop("\nD must be an f (components) by f (components) diagonal matrix.")
        if (NCOL(D) != components) stop("\nD must be an f (components) by f (components) diagonal matrix.")
        if (NROW(E) == NCOL(E)) {
            C <- E %*% D %*% t(E)
        } else {
            C <- NULL
        }
    }

    # Z (n by f)
    if (!is.null(Z)) {
        if (!is.matrix(Z)) stop("\nwhitening signal matrix Z must be an n (data rows) by f (components) matrix")
        if (NROW(Z) != n) stop("\nwhitening signal matrix Z must be an n (data rows) by f (components) matrix")
        if (NCOL(Z) != components) stop("\nwhitening signal matrix Z must be an n (data rows) by f (components) matrix")
    } else {
        if (!is.null(E)) {
            K <- solve(sqrt(D)) %*% t(E)
            Z <- M %*% t(K)
            L <- E %*% sqrt(D)
        }
    }

    if (!is.null(L)) {
        if (is.null(K)) stop("\nde-whitening matrix L cannot be provided without whitening matrix K.")
        if (NROW(L) != m) stop("\nde-whitening matrix L must be an m (signals) by f (components) matrix.")
        if (NCOL(L) != components) stop("\nde-whitening matrix L must be an m (signals) by f (components) matrix.")
    }

    # whitening matrix K (m by f)
    if (!is.null(K)) {
        if (!is.matrix(K)) stop("\nwhitening matrix K must be an f (components) by m (series) matrix.")
        if (NROW(K) != components) stop("\nwhitening matrix K must be an f (components) by m (series) matrix.")
        if (NCOL(K) != m) stop("\nwhitening matrix K must be an f (components) by m (series) matrix.")
        if (is.null(L)) {
            if (NCOL(K) == NROW(K)) {
                L <- t(solve(K))
            } else{
                stop("\ndewhitening matrix L needs to be supplied for user supplied non-symmetric K.")
            }
        }
        if (is.null(Z)) {
            Z <- M %*% t(K)
        }
    }
    if (!augment) replications <- 1
    spacings <- floor(sqrt(n))

    pca_calc_conditions <- 0
    if (!is.null(E)) pca_calc_conditions <- pca_calc_conditions + 1
    if (!is.null(D)) pca_calc_conditions <- pca_calc_conditions + 1
    whitening_calc_conditions <- 0
    if (!is.null(Z)) whitening_calc_conditions <- whitening_calc_conditions + 1
    if (!is.null(K)) whitening_calc_conditions <- whitening_calc_conditions + 1
    if (!is.null(L)) whitening_calc_conditions <- whitening_calc_conditions + 1
    if (trace) {
        cat(paste("\nsignals: ", m, sep = ""))
        cat(paste("\nsamples: ", n, "\n", sep = ""))
    }
    if (whitening_calc_conditions == 3) {
        if (trace) cat("\nPCA inputs supplied. PCA calculations not needed.\n")
    } else {
        if (pca_calc_conditions == 2) {
            if (trace) cat("\nvalues for PCA matrices supplied. PCA calculations not needed.\n")
        } else{
            tmp = .pca(M, first_eigen, last_eigen, pca_cov, trace, ...)
            E <- tmp$E
            D <- tmp$D
            C <- tmp$C
        }
        tmp <- .whitening_fun(M, E, D, trace)
        Z <- tmp$Z
        K <- tmp$K
        L <- tmp$L
    }

    sweeps <- m - 1
    if (trace) cat("\nrunning RADICAL...\n")
    solution <- .radical_recursion(k = k, sigma = sigma, samples = n, replications = replications,
                                   whitening_signal = t(Z), whitening_matrix = K, dewhitening_matrix = L,
                                   mixed_signal = t(M), mixed_mean = mu, trace = trace)
    W <- solution$W
    S <- solution$S
    U <- solution$U
    A <- solution$A
    toc = Sys.time() - tic
    return(list(A = A, W = W, S = S, U = U,
                K = K, L = L, C = C,
                Z = Z,
                D = D,
                E = E,
                mu = mu, elapsed = toc))
}

filter_ica <- function(X, demean = FALSE, W, mu)
{
    if (demean) {
        M <- sweep(X, 2, mu, FUN = "-")
    } else {
        M <- X
        mu <- rep(0, ncol(X))
    }
    # check
    S <- M %*% W + matrix(1, nrow = nrow(X)) %*% matrix(mu, nrow = 1) %*% W
    return(S)
}

#' The Fast ICA (FASTICA)
#' @description
#' The fast fixed point algorithm for independent component analysis and
#' projection pursuit based on the direct translation to R of the FastICA
#' program of the original authors at the Helsinki University of Technology.
#' @param X the data matrix (n x m) where n is the number of samples and m is the number of signals.
#' @param components the number of independent components to extract.
#' @param fun the non-linear function to use in the fixed-point algorithm.
#' @param tune the non-linear function to use for fine tuning.
#' @param method the method to use, with \dQuote{symmetric} estimating all the
#' independent component in parallel, whereas \dQuote{deflation} estimates
#' independent components one-by-one similar to projection pursuit.
#' @param tanh_par control parameter for the \dQuote{tanh} function.
#' @param gauss_par control parameter for the \dQuote{gauss} function.
#' @param step_size if this is anything other than 1, the program will use the
#' stabilized version of the algorithm (see stabilization argument).
#' @param stabilization whether to use the stabilized version of the algorithm, in
#' which case the step_size can momentarily be halved if the algorithm is stuck
#' between two points (called a stroke). Additionally if there no convergence
#' before half of the maximum number of iterations has been reached then step_size
#' will be halved for the rest of the rounds.
#' @param tol the stopping tolerance.
#' @param maxiter the maximum number of iterations.
#' @param maxiter_fine the maximum number of iterations for the fine tuning algorithm in
#' the deflation method.
#' @param sampling_ratio percent (0-1) of random samples used in one iteration
#' @param demean whether to demean the data.
#' @param pca_cov the method used to calculate the covariance matrix for PCA. Options are
#' \dQuote{ML} (maximum likelihood), \dQuote{LW} (Ledoit-Wolf) shrinkage method and
#' \dQuote{EWMA} (exponentially weighted moving average).
#' @param first_eigen This and \code{last_eigen} specify the range for eigenvalues
#' that are retained, \code{first_eigen} is the index of largest eigenvalue to be
#' retained. Making use of this option overwrites \code{components}.
#' @param last_eigen the index of the last (smallest) eigenvalue to be retained
#' and overwrites \code{component} argument.
#' @param E (optional) Eigen vector from the PCA decomposition. If this
#' is provided then \code{D} must also be provided.
#' @param D (optional) Eigen values from the PCA decomposition.
#' @param Z (optional) provided whitened signal matrix. If this is provided then
#' both \code{K} and \code{L} must also be provided.
#' @param K (optional) whitening matrix.
#' @param L (optional) de-whitening matrix.
#' @param A_init optional initialization for the mixing matrix.
#' @param seed the random seed for the random number generator.
#' @param trace whether to print out progress information.
#' @param ... additional arguments to be passed to the covariance matrix calculation.
#' For arguments passed to the \dQuote{EWMA} method, it optionally takes an additional argument
#' \code{lambda} which is the exponential decay parameter (default is 0.96).
#' The \dQuote{LW} takes an additional argument \code{shrink} which is the
#' shrinkage parameter (default is to calculate this).
#' @details
#' Steps to the general algorithm are as follows:
#' \enumerate{
#' \item Demean the data if required: \eqn{M = X - \mu}
#' \item Calculate the covariance matrix \eqn{\Sigma}  using one of the
#' methods provided.
#' \item Use an eigen decomposition to calculate the eigenvalues and eigenvectors
#' of the covariance matrix: \eqn{\Sigma = E D E'}
#' \item Select the range of eigenvalues to retain (dimensionality reduction).
#' \item Calculate the whitening matrix \eqn{K = D^{-1/2}E'} and the dewhitening
#' matrix \eqn{L = E D^{1/2}}.
#' \item Whiten the data: \eqn{Z = M K'}. Unwhitening is done by \eqn{M = Z L'}.
#' \item Run the RADICAL algorithm to calculate the rotation matrix \eqn{U},
#' the mixing matrix: \eqn{A = U^{-1} L} and the unmixing matrix \eqn{W = K' U}.
#' \item Calculate the independent components: \eqn{S = M W + \bold{1}\mu W } where
#' \bold{1} is a matrix of ones with dimension (samples x 1).
#' }
#' Notice that in calculating the mixing (A) and unmixing (W) matrices we have combined
#' the whitening (K) and un-whitening (L) matrices with the rotation matrix \eqn{U}.
#' @references
#' \insertRef{Hyvarinen1999}{tsmarch}
#' @returns A list with the following components:
#' \item{A}{the mixing matrix}
#' \item{W}{the unmixing matrix}
#' \item{S}{the independent components}
#' \item{U}{the rotation matrix}
#' \item{K}{the whitening matrix}
#' \item{L}{the dewhitening matrix}
#' \item{C}{the covariance matrix}
#' \item{Z}{the whitened signal}
#' \item{mu}{the mean of the mixed signal (X)}
#' \item{elapsed}{the time taken to run the algorithm}
#' @export
#'
fastica <- function(X, components = NCOL(X), demean = TRUE, method = "symmetric",
                    pca_cov = "ML", first_eigen = NULL, last_eigen = NULL, E = NULL, D = NULL,
                    Z = NULL, K = NULL, L = NULL, A_init = NULL,
                    fun = "tanh", tune = "none", tanh_par = 1, gauss_par = 1,
                    step_size = 1, stabilization = FALSE, tol = 1e-4, maxiter = 1000,
                    maxiter_fine = 100, sampling_ratio = 1, seed = NULL,
                    trace = FALSE, ...)
{
    tic <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    m <- NCOL(X)
    n <- NROW(X)
    if (m > n) warning("\nsignal matrix (X) may be orientated in the wrong way.")
    pca_cov <- match.arg(pca_cov[1], c("ML", "LW", "EWMA"))
    method <- match.arg(method[1], c("symmetric","deflation"))
    fun <- match.arg(fun[1], c("pow3", "tanh", "gauss", "skew"))
    stabilization <- as.logical(stabilization)
    step_size <- as.numeric(step_size)
    tanh_par <- as.numeric(tanh_par)
    gauss_par <- as.numeric(gauss_par)
    maxiter <- as.integer(abs(maxiter))
    maxiter_fine <- as.integer(abs(maxiter_fine))
    tol <- as.numeric(abs(tol))
    trace <- as.logical(trace)
    if (demean) {
        mu <- apply(X, 2, mean)
        M <- sweep(X, 2, mu, FUN = "-")
    } else {
        mu <- rep(0, m)
        M <- X
    }
    if (components > m) stop("\ncannot choose more components than signals.")
    if (is.null(first_eigen)) first_eigen <- 1
    if (is.null(last_eigen)) last_eigen <- components
    if (!is.null(first_eigen)) {
        if (is.null(last_eigen)) {
            last_eigen <- last_eigen + components
            if (last_eigen > m) stop("\nfirst_eigen + components is greater than number of signals.")
        } else {
            if (last_eigen > m) stop("\nlast_eigen greater than number of signals.")
            components <- last_eigen - first_eigen + 1
            if (first_eigen < 1) stop("\ncannot have zero factors in PCA.")
        }
    }
    # PCA E (m x f)
    if (!is.null(E)) {
        if (is.null(D)) stop("\nD required when passing E.")
        if (!is.matrix(E)) stop("\nE must be an m (signals) x f (components) matrix.")
        if (NROW(E) != m) stop("\nE must be an m (signals) x f (components) matrix.")
        if (NCOL(E) != components) stop("\nE must be an m (signals) x f (components) matrix.")
    }

    # PCA D (f by f)
    if (!is.null(D)) {
        if (is.null(D)) stop("\nE required when passing D.")
        if (!is.matrix(D)) stop("\nD must be an f (components) by f (components) diagonal matrix.")
        if (NROW(D) != components) stop("\nD must be an f (components) by f (components) diagonal matrix.")
        if (NCOL(D) != components) stop("\nD must be an f (components) by f (components) diagonal matrix.")
        if (NROW(E) == NCOL(E)) {
            C <- E %*% D %*% t(E)
        } else {
            C <- NULL
        }
    }
    # Z (n by f)
    if (!is.null(Z)) {
        if (!is.matrix(Z)) stop("\nwhitening signal matrix Z must be an n (data rows) by f (components) matrix")
        if (NROW(Z) != n) stop("\nwhitening signal matrix Z must be an n (data rows) by f (components) matrix")
        if (NCOL(Z) != components) stop("\nwhitening signal matrix Z must be an n (data rows) by f (components) matrix")
    } else {
        if (!is.null(E)) {
            K <- solve(sqrt(D)) %*% t(E)
            Z <- M %*% t(K)
            L <- E %*% sqrt(D)
        }
    }
    if (!is.null(L)) {
        if (is.null(K)) stop("\nde-whitening matrix L cannot be provided without whitening matrix K.")
        if (NROW(L) != m) stop("\nde-whitening matrix L must be an m (signals) by f (components) matrix.")
        if (NCOL(L) != components) stop("\nde-whitening matrix L must be an m (signals) by f (components) matrix.")
    }
    # whitening matrix K (m by f)
    if (!is.null(K)) {
        if (!is.matrix(K)) stop("\nwhitening matrix K must be an f (components) by m (series) matrix.")
        if (NROW(K) != components) stop("\nwhitening matrix K must be an f (components) by m (series) matrix.")
        if (NCOL(K) != m) stop("\nwhitening matrix K must be an f (components) by m (series) matrix.")
        if (is.null(L)) {
            if (NCOL(K) == NROW(K)) {
                L <- t(solve(K))
            } else{
                stop("\ndewhitening matrix L needs to be supplied for user supplied non-symmetric K.")
            }
        }
        if (is.null(Z)) {
            Z <- M %*% t(K)
        }
    }

    pca_calc_conditions <- 0
    if (!is.null(E)) pca_calc_conditions <- pca_calc_conditions + 1
    if (!is.null(D)) pca_calc_conditions <- pca_calc_conditions + 1
    whitening_calc_conditions <- 0
    if (!is.null(Z)) whitening_calc_conditions <- whitening_calc_conditions + 1
    if (!is.null(K)) whitening_calc_conditions <- whitening_calc_conditions + 1
    if (!is.null(L)) whitening_calc_conditions <- whitening_calc_conditions + 1
    if (trace) {
        cat(paste("\nsignals: ", m, sep = ""))
        cat(paste("\nsamples: ", n, "\n", sep = ""))
    }
    if (whitening_calc_conditions == 3) {
        if (trace) cat("\nPCA inputs supplied. PCA calculations not needed.\n")
    } else {
        if (pca_calc_conditions == 2) {
            if (trace) cat("\nvalues for PCA matrices supplied. PCA calculations not needed.\n")
        } else{
            tmp = .pca(M, first_eigen, last_eigen, pca_cov, trace, ...)
            E <- tmp$E
            D <- tmp$D
            C <- tmp$C
        }
        tmp <- .whitening_fun(M, E, D, trace)
        Z <- tmp$Z
        K <- tmp$K
        L <- tmp$L
    }

    sampling_ratio <- abs(sampling_ratio)
    if (sampling_ratio > 1) {
        sampling_ratio <- 1
        warning("\nsetting sampling_ratio to 1")
    } else {
        if ((sampling_ratio * n) < 1000) {
            sampling_ratio <- min(1000/n, 1)
        }
    }
    if (!is.null(tune)) {
        tune <- match.arg(tune[1], c("none","pow3", "tanh", "gauss", "skew"))
        if (tune == "none") enable_tuning <- FALSE else enable_tuning <- TRUE
    } else {
        tune <- "none"
        enable_tuning <- FALSE
    }
    # Checking the value for nonlinearity.
    fun <- switch(tolower(fun),
                  "pow3"  = 10,
                  "tanh"  = 20,
                  "gauss" = 30,
                  "skew"  = 40)
    if (sampling_ratio != 1) fun <- fun + 2
    if (step_size  != 1) fun <- fun + 1

    tune <- switch(tolower(tune),
                   "pow3"  = 10 + 1,
                   "tanh"  = 20 + 1,
                   "gauss" = 30 + 1,
                   "skew"  = 40 + 1,
                   "none"  = if (step_size != 1) fun else fun + 1)

    if (stabilization) {
        stabilization <- TRUE
    } else {
        if (step_size != 1) stabilization <- TRUE else stabilization <- FALSE
    }

    if (is.null(A_init)) {
        init_state_mode <- 0
        if (trace) cat("\nusing random initialization for A.")
    } else {
        if (NROW(A_init) != NCOL(K)) {
            init_state_mode <- 0
            A_init <- NULL
            warning("\nsize of initial guess is incorrect. Using random initial guess.")
        } else {
            init_state_mode <- 1
            if (NCOL(A_init) < components) {
                warning(paste("\ninitial guess only for first ", NCOL(A_init)," components. Using random initial A_init for others.", sep = ""))
                set.seed(seed)
                A_init[,NCOL(A_init) + (1:components)] <- matrix(runif(m * (components - NCOL(A_init))) - 0.5, ncol = components - NCOL(A_init))
            } else{
                if (NCOL(A_init) > components) warning(paste("\nInitial guess too large. The excess column are dropped.", sep = ""))
                A_init <- A_init[, 1:components]
            }
            if (trace) cat("\nusing A_init.")
        }
    }
    if (trace) cat("\nstarting ICA calculation...\n")
    solution <- switch(method,
                       "symmetric" = fp_symmetric(Z, fun, tune, enable_tuning, init_state_mode, K, L, A_init, maxiter, tol, stabilization, step_size, sampling_ratio, tanh_par, gauss_par, trace),
                       "deflation" = fp_deflation(Z, fun, tune, enable_tuning, init_state_mode, K, L, A_init, maxiter, maxiter_fine,
                                                  tol, stabilization, step_size, sampling_ratio, tanh_par, gauss_par, trace))
    W <- t(solution$W)
    U <- solution$U
    A <- t(solution$A)
    L <- t(L)
    S <- M %*% W + matrix(1, nrow = n) %*% (mu %*% W)
    toc = Sys.time() - tic
    return(list(A = A, W = W, S = S, U = U,
                K = K, L = L, C = C,
                Z = Z,
                D = D,
                E = E,
                mu = mu, elapsed = toc))
}

fp_symmetric <- function(X, fun, tune, enable_tuning, init_state_mode = 0, K, L, A_init,
                          maxiter, tol, stabilization, step_size, sampling_ratio,
                          tanh_par, gauss_par, trace)
{
    components <- NCOL(X)
    n <- NROW(X)
    initial_step_size <- step_size
    step_size_k <- 0.01
    stroke <- 0
    not_fine <- TRUE
    long <- FALSE
    n_fun <- fun
    A <- matrix(0, components, components)
    if (init_state_mode == 0) {
        # Take random orthonormal initial vectors.
        U <- .orthogonal(matrix(rnorm(components * components), ncol = components, byrow = TRUE))
    } else {
        # Use the given initial vector as the initial state
        U <- K %*% A_init
    }
    U_old_2 <- U_old <- matrix(0, NROW(U), NCOL(U))
    # This is the actual fixed-point iteration loop.
    for (i in 1:(maxiter + 1)) {
        if (i == maxiter + 1) {
            if (trace) cat(paste("No convergence after ", maxiter," steps\n", sep = ""))
            if (!is.null(U)) {
                U <- U * .sqrtm(Re(solve(t(U) %*% U)))
                W <- t(U) %*% K
                A <- L %*% U
            } else {
                W <- NULL
                A <- NULL
            }
            # exit condition
            return(list(A = A, W = W))
        }
        # Symmetric orthogonalization.
        U <- U %*% Re(.sqrtm(solve(t(U) %*% U)))
        # Test for termination condition. Note that we consider opposite
        # directions here as well.
        min_abs_cos <- min(abs(diag(t(U) %*% U_old)))
        min_abs_cos_2 <- min(abs(diag(t(U) %*% U_old_2)))
        if ((1 - min_abs_cos) < tol) {
            if (enable_tuning && not_fine) {
                if (trace) cat("\ninitial convergence, fine-tuning: \n")
                not_fine <- FALSE
                n_fun <- tune
                step_size <- step_size_k * initial_step_size
                U_old <- U_old_2 <- matrix(0, NROW(U), NCOL(U))
            } else {
                if (trace) cat(paste("\nconvergence after ", i, "steps\n", sep = ""))
                # Calculate the de-whitened vectors.
                A <- L %*% U
                break()
            }
        }
        if (stabilization) {
            if (stroke == 0.0) {
                if ((1 - min_abs_cos_2) < tol) {
                    if (trace) cat("\nstroke!\n")
                    stroke <- step_size
                    step_size <- 0.5 * step_size
                    if (n_fun %% 2 == 0) n_fun <- n_fun + 1
                }
            } else {
                step_size <- stroke
                stroke <- 0.0
                if ((step_size == 1) && (n_fun %% 2 != 0)) n_fun <- n_fun - 1
            }
            if (!long && (i > maxiter/2)) {
                if (trace) cat("\ntaking long (reducing step size)\n")
                long <- TRUE
                step_size <- 0.5 * step_size
                if (n_fun %% 2 == 0) n_fun <- n_fun + 1
            }
        }
        U_old_2 <- U_old
        U_old <- U
        if (trace) {
            if (i == 1) {
                cat(paste("\nstep no. ", i,"\n", sep = ""))
            } else {
                cat(paste("\nstep no. ", i,", change in value of estimate ", round(1 - min_abs_cos,3),"\n", sep = ""))
            }
        }
        U = .fs_method(func_no = n_fun, X = t(X), B = U, step_size = step_size, n = n, sampling_ratio = sampling_ratio,
                       tanh_par = tanh_par, gauss_par = gauss_par)
    }
    W <- t(U) %*% K
    return(list(W = W, A = A, U = U))
}


fp_deflation <- function(X, fun, tune, enable_tuning, init_state_mode = 0, K, L, A_init,
                         maxiter, maxiter_fine, tol, stabilization, step_size, sampling_ratio,
                         tanh_par, gauss_par, trace)
{
    components <- NCOL(X)
    n <- NROW(X)
    step_size_k <- 0.01
    initial_step_size <- step_size
    B <- matrix(0, components, components)
    A <- matrix(0, NROW(L), NCOL(L))
    W <- matrix(0, NCOL(L), NROW(L))
    num_failures <- 0
    failure_limit <- 5
    j <- 1
    while (j <= components) {
        step_size <- initial_step_size
        stroke <- 0
        not_fine <- TRUE
        long <- FALSE
        n_fun <- fun
        end_fine_tuning <- 0
        if (trace) cat(paste0("\nIC_", j,"..."))
        if (is.null(A_init)) {
            w <- matrix(rnorm(components), ncol = 1)
        } else{
            w <- K %*% A_init[,j]
        }
        w <- w - B %*% t(B) %*% w
        w <- w / .fnorm(w)
        w_old_2 <- w_old <- .zeros(NROW(w), NCOL(w))
        i <- 1
        gabba <- 1
        while (i <= (maxiter + gabba)) {
            w <- w - B %*% t(B) %*% w
            w <- w / .fnorm(w)
            if (not_fine) {
                if (i == (maxiter + 1)) {
                    if (trace) cat(paste0("\ncomponent number ", j, " did not converge in ", maxiter, " iterations."))
                    j <- j - 1
                    num_failures <- num_failures + 1
                    if (num_failures > failure_limit) {
                        if (trace)  cat(paste0("\ntoo many failures to converge (",  num_failures, "). Giving up."))
                        if (j == 0) {
                            A <- NULL
                            W <- NULL
                        }
                        return(list(W = W, A = A))
                    }
                }
            } else{
                if (i >= end_fine_tuning) w_old <- w
            }
            if ((.fnorm(w - w_old) < tol) | (.fnorm(w + w_old) < tol)) {
                if (enable_tuning && not_fine) {
                    if (trace) cat("\ninitial convergence, fine-tuning")
                    not_fine <- FALSE
                    gabba <- maxiter_fine
                    w_old_2 <- w_old <- .zeros(NROW(w), NCOL(w))
                    n_fun <- tune
                    step_size <- step_size_k * initial_step_size
                    end_fine_tuning <- maxiter_fine + i
                } else {
                    num_failures <- 0
                    # Save the vector
                    B[,j] <- w
                    # Calculate the de-whitened vector.
                    A[,j] <- L %*% w
                    # Calculate ICA filter.
                    W[j,] <- t(w) %*% K
                    if (trace) cat(paste0("\ncomputed (", i, " steps)"))
                    break()
                }
            }
            if (stabilization) {
                if (stroke == 0 && ((.fnorm(w - w_old_2) < tol) | (.fnorm(w + w_old_2) < tol))) {
                    stroke <- step_size
                    if (trace) cat("\nstroke!")
                    step_size <- 0.5 * step_size
                    if ((n_fun %% 2) == 0) n_fun <- n_fun + 1
                }
                if (stroke > 0) {
                    step_size <- stroke
                    stroke <- 0
                    if ((step_size == 1) && ((n_fun %% 2) != 0)) n_fun <- n_fun - 1
                }
                if (not_fine && !long && (i > maxiter/2)) {
                    if (trace) cat("\ntaking too long (reducing step size")
                    long <- TRUE
                    step_size <- 0.5 * step_size
                    if ((n_fun %% 2) == 0) n_fun <- n_fun + 1
                }
            }
            w_old_2 <- w_old
            w_old <- w
            w <- .fd_method(n_fun, t(X), w, step_size, n, sampling_ratio, tanh_par, gauss_par)
            w <- w / .fnorm(w)
            i <- i + 1
        }
        j = j + 1
    }
    U <- t(A) %*% t(K)
    return(list(W = W, A = A, U = t(U)))
}
