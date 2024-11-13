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
        mu <- colMeans(X)
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

filter_radical <- function(X, demean = FALSE, W, mu)
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

.pca <- function(X, first_eigen = 1, last_eigen = dim(X)[1], pca_cov, trace, ...)
{
    old_dimension <- dim(X)[2]
    Y <- X
    covariance_matrix <- switch(pca_cov,
                               "ML" = (t(Y) %*% Y)/dim(Y)[1],
                               "LW" = lw_covariance(Y, demean = FALSE,...),
                               "EWMA" = ewma_covariance(Y, demean = FALSE, ...))
    ed <- eigen(covariance_matrix)
    D <- diag(ed$values)
    E <- ed$vectors
    rank_tol <- 1e-7
    max_last_eigen <- sum(diag(D) > rank_tol)
    if (max_last_eigen == 0) {
        stop("\nEigenvalues of the calculated covariance matrix are all smaller than tolerance of 1e-7. Try rescaling the data matrix.")
    }
    eigen_values <- sort(diag(D), decreasing  = TRUE)
    if (last_eigen > max_last_eigen) {
        last_eigen <- max_last_eigen
        if (trace) cat(paste("\nDimension reduced to ",last_eigen - first_eigen + 1," due to the singularity of covariance matrix\n", sep = ""))
    } else {
        if (trace) {
            if (old_dimension == (last_eigen - first_eigen + 1)) {
                cat("\ndimension not reduced.\n")
            } else {
                cat("\nreducing dimension...\n")
            }
        }
    }
    # Drop the smaller eigenvalues
    if (last_eigen < old_dimension) {
        lower_limit_value <- (eigen_values[last_eigen] + eigen_values[last_eigen + 1]) / 2
    } else{
        lower_limit_value <- eigen_values[old_dimension] - 1
    }
    lower_cols <- diag(D) > lower_limit_value
    # Drop the larger eigenvalues
    if (first_eigen > 1) {
        higher_limit_value = (eigen_values[first_eigen - 1] + eigen_values[first_eigen]) / 2
    } else {
        higher_limit_value = eigen_values[1] + 1
    }
    higher_cols <- diag(D) < higher_limit_value

    # Combine the results from above
    selected_cols <- lower_cols & higher_cols

    # print some info for the user
    if (trace) cat(paste("selected ", sum(selected_cols), "dimensions.\n", sep = ""))

    if (sum(selected_cols) != (last_eigen - first_eigen + 1)) {
        stop("\nselected wrong number of dimensions.\n")
    }

    if (trace) {
        cat(paste("smallest remaining (non-zero) eigenvalue ", round(eigen_values[last_eigen], 5), "\n", sep = ""))
        cat(paste("largest remaining (non-zero) eigenvalue ", round(eigen_values[first_eigen], 5), "\n", sep = ""))
        cat(paste("sum of removed eigenvalues ", round(sum(diag(D) * (!selected_cols)), 5), "\n", sep = ""))
    }

    # Select the colums which correspond to the desired range of eigenvalues.
    E <- .sel_col(E, selected_cols)
    D <- .sel_col(t(.sel_col(D, selected_cols)), selected_cols)
    if (trace) {
        sum_all <- sum(eigen_values)
        sum_used <- sum(diag(D))
        sum_retained <- (sum_used / sum_all) * 100
        cat(paste(round(sum_retained, 2), " % of (non-zero) eigenvalues retained.\n", sep = ""))
    }
    return(list(E = E, D = D, C = covariance_matrix))
}

.sel_col <- function(x, mask)
{
    # Selects the columns of the matrix that marked by one in the given vector.
    # The mask is a column vector.
    use <- numeric()
    if (length(mask) != NCOL(x)) {
        stop("\nthe mask vector must be of length NCOL(x).\n")
    }
    idx <- 0
    for (i in 1:length(mask)) {
        if (mask[i] == 1) {
            use[idx + 1] <- i
            idx <- idx + 1
        }
    }
    out = x[, use, drop = FALSE]
    return(out)
}

.whitening_fun <- function(x, E, D, trace)
{
    if (any(diag(D) < 0)) stop("\nnegative eigenvalues computed from the covariance matrix")
    K <- solve(sqrt(D)) %*% t(E)
    L <- E %*% sqrt(D)
    if (trace) cat("Whitening...\n")
    Z <- x %*% t(K)
    if (any(is.complex(Z))) stop("\nwhitened matrix has imaginary values.")
    return(list(Z = Z, K = K, L = L))
}

lw_covariance <- function(X, shrink = -1, demean = FALSE, trace) {
    n <- NROW(X)
    m <- NCOL(X)
    mu <- colMeans(X)
    if (demean) X <- sweep(X, 2, mu, FUN = "-")
    # compute sample covariance matrix
    sample_covariance <- (t(X) %*% X)/n
    # compute prior
    mean_var <- mean(diag(sample_covariance))
    prior <- mean_var * diag(1, m, m)
    if (shrink == -1) {
        # compute shrinkage parameters
        # p in paper
        Y <- X^2
        phi_mat <- (t(Y) %*% Y)/n - 2 * (t(X) %*% X) * sample_covariance/n + sample_covariance^2
        phi <- sum(apply(phi_mat, 1, "sum"))
        # c in paper
        cgamma <- norm(sample_covariance - prior, 'F')^2
        # shrinkage constant
        kappa <- phi/cgamma
        shrinkage <- max(0, min(1, kappa/n))
        if (trace) cat(paste("shrinkage parameter: ", shrinkage, "\n", sep = ""))
    } else {
        shrinkage <- shrink
    }
    sigma <- shrinkage * prior + (1 - shrinkage) * sample_covariance
    return(sigma)
}

ewma_covariance <- function(X, lambda = 0.96, demean = FALSE)
{
    n <- NROW(X)
    i <- (0:(n - 1))
    ewma_wt <- lambda^i
    ewma_wt <- ewma_wt/sum(ewma_wt)
    covariance <- cov.wt(X, wt = rev(ewma_wt), center = demean)$cov
    return(covariance)
}
