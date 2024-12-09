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

.orthogonal <- function(x)
{
    tmp <- svd(x)
    u <- tmp$u
    s <- tmp$d
    m <- dim(x)[1]
    n <- dim(x)[2]
    if (n == 1 || m == 1)  s <- s[1]
    tol <- max(m, n) * s[1] * .Machine$double.eps
    rnk <- sum(s > tol)
    ans <- u[,1:rnk]
    return(ans)
}

.sqrtm <- function(x)
{
    tmp <- svd(x)
    sqrt_x <- tmp$u %*% sqrt(diag(tmp$d)) %*% t(tmp$u)
    return(sqrt_x)
}

# frobenius norm
.fnorm = function(x)
{
    sqrt(sum(as.numeric(x)^2))
}

.zeros <- function(n = 1, m = 1)
{
    if (missing(m)) m <- n
    return(matrix(0, nrow = n, ncol = m))
}

.ones = function(n = 1, m = 1)
{
    if (missing(m)) m <- n
    return(matrix(1, nrow = n, ncol = m))
}

.samples <- function(max_n, p)
{
    ans <- numeric()
    while (length(ans) == 0) {
        ans <- which(runif(max_n) < p)
    }
    return(ans)
}

.fs_method <- function(func_no, X, B, step_size, n, sampling_ratio, tanh_par, gauss_par)
{
    ans <- switch(as.character(func_no),
                  "10" = .fs_10(X, B,  step_size, n),
                  "11" = .fs_11(X, B,  step_size, n),
                  "12" = .fs_12(X, B,  step_size, n,  sampling_ratio),
                  "13" = .fs_13(X, B,  step_size, n,  sampling_ratio),
                  "20" = .fs_20(X, B,  step_size, tanh_par,  n),
                  "21" = .fs_21(X, B,  step_size, tanh_par),
                  "22" = .fs_22(X, B,  step_size, tanh_par,  n, sampling_ratio),
                  "23" = .fs_23(X, B,  step_size, tanh_par,  n, sampling_ratio),
                  "30" = .fs_30(X, B,  step_size, gauss_par, n, sampling_ratio),
                  "31" = .fs_31(X, B,  step_size, gauss_par),
                  "32" = .fs_32(X, B,  step_size, gauss_par, n, sampling_ratio),
                  "33" = .fs_33(X, B,  step_size, gauss_par, n, sampling_ratio),
                  "40" = .fs_40(X, B,  step_size, n),
                  "41" = .fs_41(X, B,  step_size),
                  "42" = .fs_42(X, B,  step_size, n,  sampling_ratio),
                  "43" = .fs_43(X, B,  step_size, n,  sampling_ratio))
    return(ans)
}

.fd_method <- function(func_no, X, B, step_size, n, sampling_ratio, tanh_par, gauss_par)
{
    ans <- switch(as.character(func_no),
                  "10" = .fd_10(X, B,  step_size, n),
                  "11" = .fd_11(X, B,  step_size, n),
                  "12" = .fd_12(X, B,  step_size, n,  sampling_ratio),
                  "13" = .fd_13(X, B,  step_size, n,  sampling_ratio),
                  "20" = .fd_20(X, B,  step_size, tanh_par,  n),
                  "21" = .fd_21(X, B,  step_size, tanh_par),
                  "22" = .fd_22(X, B,  step_size, tanh_par,  n, sampling_ratio),
                  "23" = .fd_23(X, B,  step_size, tanh_par,  n, sampling_ratio),
                  "30" = .fd_30(X, B,  step_size, gauss_par, n),
                  "31" = .fd_31(X, B,  step_size, gauss_par),
                  "32" = .fd_32(X, B,  step_size, gauss_par, n, sampling_ratio),
                  "33" = .fd_33(X, B,  step_size, gauss_par, n, sampling_ratio),
                  "40" = .fd_40(X, B,  step_size, n),
                  "41" = .fd_41(X, B,  step_size, n),
                  "42" = .fd_42(X, B,  step_size, n,  sampling_ratio),
                  "43" = .fd_43(X, B,  step_size, n,  sampling_ratio))
    return(ans)
}

# pow3
.fs_10 <- function(X, B, step_size, n)
{
    ans <- (X %*% ((t(X) %*% B)^3)) / n - 3.0 * B
    return(ans)
}

.fs_11 <- function(X, B, step_size, n)
{
    Y <- t(X) %*% B
    G <- Y^3
    Beta <- as.numeric(apply(Y * G, 2, "sum"))
    D <- diag(1 / (Beta - 3 * n))
    ans <- B + step_size * B %*% (t(Y) %*% G - diag(Beta)) %*% D
    return(ans)
}

.fs_12 <- function(X, B,step_size, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    ans <- (X_subset %*% (( t(X_subset) %*% B)^3)) / NCOL(X_subset) - 3 * B
    return(ans)
}

.fs_13 <- function(X, B,step_size, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- t(X_subset) %*% B
    G <- Y^3
    Beta <- as.numeric(apply(Y %*% G, 2, "sum"))
    D <- diag(1 / (Beta - 3 * NCOL(t(Y))))
    ans <- B + step_size * B %*% (t(Y) %*% G - diag(Beta)) %*% D
    return(ans)
}

# tanh
.fs_20 <- function(X, B, step_size, tanh_par, n)
{
    hyp_tan <- tanh(tanh_par * t(X) %*% B)
    ans <- X %*% hyp_tan / n - .ones(NROW(B),1) %*% apply(1 - hyp_tan^2, 2, "sum") * B / n * tanh_par
    return(ans)
}


.fs_21 <- function(X, B, step_size, tanh_par)
{
    Y <- t(X) %*% B
    hyp_tan <- tanh(tanh_par * Y)
    Beta <- apply(Y * hyp_tan, 2, "sum")
    D <- diag(1/(Beta - tanh_par * apply(1 - hyp_tan^2, 2, "sum")))
    ans <- B + step_size * B %*% (t(Y) %*% hyp_tan - diag(Beta)) %*% D
    return(ans)
}

.fs_22 <- function(X, B, step_size, tanh_par, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    hyp_tan <- tanh(tanh_par * t(X_subset) %*% B)
    ans <- X_subset %*% hyp_tan / NCOL(X_subset) - .ones(NROW(B),1) %*% apply(1 - hyp_tan^2, 2, "sum") * B / NCOL(X_subset) * tanh_par
    return(ans)
}

.fs_23 <- function(X, B, step_size, tanh_par, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- t(X_subset) %*% B
    hyp_tan <- tanh(tanh_par * Y)
    Beta <- apply(Y * hyp_tan, 2, "sum")
    D <- diag(1 / (Beta - tanh_par * apply(1 - hyp_tan^2, 2, "sum")))
    ans <- B + step_size * B %*% (t(Y) %*% hyp_tan - diag(Beta)) %*% D
    return(ans)
}

# gauss
.fs_30 <- function(X, B, step_size, gauss_par, n, sampling_ratio)
{
    U <- t(X) %*% B
    U_squared <- U^2
    ex <- exp(-gauss_par * U_squared / 2)
    gauss <- U * ex
    d_gauss <- (1 - gauss_par * U_squared) * ex
    ans <- X %*% gauss / n - .ones(NROW(B),1) %*% apply(d_gauss, 2, "sum") * B / n
    return(ans)
}

.fs_31 <- function(X, B, step_size, gauss_par)
{
    Y <- t(X) %*% B
    ex <- exp(-gauss_par * (Y^2) / 2)
    gauss <- Y * ex
    Beta <- apply(Y * gauss, 2, "sum")
    D <- diag(1 / (Beta - apply((1 - gauss_par * (Y^2)) * ex, 2, "sum")))
    ans <- B + step_size * B %*% (t(Y) %*% gauss - diag(Beta)) %*% D
    return(ans)
}


.fs_32 <- function(X, B, step_size, gauss_par, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    U <- t(X_subset) %*% B
    U_squared <- U^2
    ex <- exp(-gauss_par * U_squared / 2)
    gauss <- U * ex
    d_gauss <- (1 - gauss_par * U_squared) * ex
    ans <- X_subset %*% gauss / NCOL(X_subset) - .ones(NROW(B), 1) * apply(d_gauss, 2, "sum") * B / NCOL(X_subset)
    return(ans)
}

.fs_33 <- function(X, B, step_size, gauss_par, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- t(X_subset) %*% B
    ex <- exp(-gauss_par * (Y^2) / 2)
    gauss <- Y * ex
    Beta <- apply(Y * gauss, 2, "sum")
    D <- diag(1 / (Beta - sum((1 - gauss_par * (Y^2)) * ex)))
    ans <- B + step_size * B %*% (t(Y) %*% gauss - diag(Beta)) %*% D
    return(ans)
}

# skew
.fs_40 <- function(X, B, step_size, n)
{
    ans <- (X %*% ((t(X) %*% B)^2)) / n
    return(ans)
}

.fs_41 <- function(X, B, step_size)
{
    Y <- t(X) %*% B
    G <- Y^2
    Beta <- apply(Y * G, 2, "sum")
    D <- diag(1 / (Beta))
    ans <- B + step_size * B %*% (t(Y) %*% G - diag(Beta)) %*% D
    return(ans)
}

.fs_42 = function(X, B, step_size, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    ans <- (X_subset %*% ((t(X_subset) %*% B)^2)) / NCOL(X_subset)
    return(ans)
}

.fs_43 <- function(X, B, step_size, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- t(X_subset) %*% B
    G <- Y^2
    Beta <- apply(Y * G, 2, "sum")
    D <- diag(1 / (Beta))
    ans <- B + step_size * B %*% (t(Y) %*% G - diag(Beta)) %*% D
    return(ans)
}

.fd_10 <- function(X, B, step_size, n)
{
    ans <- (X %*% ((t(X) %*% B)^3)) / n - 3 * B
    return(ans)
}

.fd_11 <- function(X, B, step_size, n)
{
    Y <- (X %*% ((t(X) %*% B)^3)) / n
    Beta <- as.numeric(t(B) %*% Y)
    ans <- B - step_size * (Y - Beta * B) / (3 - Beta)
    return(ans)
}

.fd_12 <- function(X, B, step_size, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    ans <- (X_subset %*% ((t(X_subset) %*% B)^3)) / NCOL(X_subset) - 3 * B
    return(ans)
}

.fd_13 <- function(X, B, step_size, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- (X_subset %*% ((t(X_subset) %*% B)^3))/NCOL(X_subset)
    Beta <- as.numeric(t(B) %*% Y)
    ans <- B - step_size * (Y - Beta * B) / (3 - Beta)
    return(ans)
}

.fd_20 <- function(X, B,  step_size, tanh_par,  n)
{
    hyp_tan <- tanh(tanh_par * t(X) %*% B)
    ans <- (X %*% hyp_tan - tanh_par * sum(1 - hyp_tan^2) * B) / n
    return(ans)
}

.fd_21 <- function(X, B,  step_size, tanh_par)
{
    hyp_tan <- tanh(tanh_par * t(X) %*% B)
    Beta <- as.numeric(t(B) %*% X %*% hyp_tan)
    ans <- B - step_size * ((X %*% hyp_tan - Beta * B)/(tanh_par * sum(1 - hyp_tan^2) - Beta))
    return(ans)
}

.fd_22 <- function(X, B,  step_size, tanh_par,  n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    hyp_tan <- tanh(tanh_par * t(X_subset) %*% B)
    ans <- (X_subset %*% hyp_tan - tanh_par * sum(1 - hyp_tan^2) * B)/NCOL(X_subset)
    return(ans)
}

.fd_23 <- function(X, B,  step_size, tanh_par,  n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    hyp_tan <- tanh(tanh_par * t(X_subset) %*% B)
    Beta <- t(B) %*% X_subset %*% hyp_tan
    ans <- B - step_size * ((X_subset %*% hyp_tan - Beta * B)/(tanh_par * sum(1 - hyp_tan^2) - Beta))
    return(ans)
}

.fd_30 <- function(X, B,  step_size, gauss_par, n)
{
    Y <- t(X) %*% B
    Y_2 <- Y^2
    ex <- exp(-gauss_par * Y_2/2)
    gauss <- Y * ex
    d_gauss <- (1 - gauss_par * Y_2) * ex
    ans <- (X %*% gauss - sum(d_gauss) * B) / n
    return(ans)
}

.fd_31 <- function(X, B,  step_size, gauss_par)
{
    Y <- t(X) %*% B
    Y_2 <- Y^2
    ex <- exp(-gauss_par * Y_2/2)
    gauss <- Y * ex
    d_gauss <- (1 - gauss_par * Y_2)*ex
    Beta <- as.numeric(t(B) %*% X %*% gauss)
    ans <- B - step_size * ((X %*% gauss - Beta * B) / (sum(d_gauss) - Beta))
    return(ans)
}

.fd_32 <- function(X, B,  step_size, gauss_par, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- t(X_subset) %*% B
    Y_2 <- Y^2
    ex <- exp(-gauss_par * Y_2/2)
    gauss <- Y * ex
    d_gauss <- (1 - gauss_par * Y_2) * ex
    ans <- (X_subset %*% gauss - sum(d_gauss) * B)/NCOL(X_subset)
    return(ans)
}

.fd_33 <- function(X, B,  step_size, gauss_par, n, sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- t(X_subset) %*% B
    Y_2 <- Y^2
    ex <- exp(-gauss_par * Y_2/2)
    gauss <- Y * ex
    d_gauss <- (1 - gauss_par * Y_2) * ex
    Beta <- as.numeric(t(B) %*% X_subset %*% gauss)
    ans <- B - step_size * ((X_subset %*% gauss - Beta * B)/(sum(d_gauss) - Beta))
    return(ans)
}

.fd_40 <- function(X, B,  step_size, n)
{
    ans <- (X %*% ((t(X) %*% B)^2)) / n
    return(ans)
}

.fd_41 <- function(X, B,  step_size, n)
{
    Y <- (X %*% ((t(X) %*% B)^2))/n
    Beta <- as.numeric(t(B) %*% Y)
    ans <- B - step_size * (Y - Beta * B)/(-Beta)
    return(ans)
}

.fd_42 <- function(X, B,  step_size, n,  sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    ans <- (X_subset %*% ((t(X_subset) %*% B)^2))/NCOL(X_subset)
    return(ans)
}

.fd_43 <- function(X, B,  step_size, n,  sampling_ratio)
{
    X_subset <- X[, .samples(n, sampling_ratio)]
    Y <- (X_subset %*% ((t(X_subset) %*% B)^2))/NCOL(X_subset)
    Beta <- as.numeric(t(B) %*% Y)
    ans <- B - step_size * (Y - Beta * B)/(-Beta)
    return(ans)
}
