.estimate_gogarch_components <- function(components, object) {
    components <- xts(components, object$target$index)
    garch_model <- object$garch$model
    garch_order <- object$garch$order
    lambda_range <- object$garch$lambda_range
    shape_range <- object$garch$shape_range
    distribution <- object$distribution
    rmod <- future_lapply(1:ncol(components), function(i){
        lower <- upper <- parameter <- NULL
        spec <- garch_modelspec(components[,i], model = garch_model, constant = FALSE, order = garch_order,
                                distribution = distribution)
        if (distribution == "gh") {
            spec$parmatrix[parameter == "lambda", lower := lambda_range[1]]
            spec$parmatrix[parameter == "lambda", upper := lambda_range[2]]
            spec$parmatrix[parameter == "shape", lower := shape_range[1]]
            spec$parmatrix[parameter == "shape", upper := shape_range[2]]
        }
        mod <- estimate(spec, keep_tmb = TRUE)
        return(mod)
    }, future.packages = c("xts","tsgarch"), future.seed = TRUE)
    rmod <- eval(rmod)
    rmod <- to_multi_estimate(rmod)
    names(rmod) <- paste0("ica_component.",1:ncol(components))
    return(rmod)
}


.gogarch_log_likelihood <- function(garch_models, K)
{
    L <- sapply(garch_models, function(x) as.numeric(logLik(x)))
    if (is_square(K)) {
        out <- sum(L) + log(abs(det(K)))
    } else {
        out <- sum(L) + log(abs(det(K %*% t(K))))
    }
    return(out)
}

.coskewness_pair_index <- function(m, k)
{
    idx <- 1:m
    idx <- idx  + ((k[3] - 1) * m)
    idx <- idx[1] + (k[2] - 1)
    return(idx)
}

.cokurtosis_pair_index <- function(m, k)
{
    idx <- 1:m
    idx <- idx  + ((k[4] - 1) * m * m)
    idx <- idx[1] + ((k[3] - 1) * m)
    idx <- idx[1] + (k[2] - 1)
    return(idx)
}


fold3d <- function(unfolded, p) {
    dims <- dim(unfolded)
    k <- dims[length(dims)]
    n <- dims[1]
    # Determine the dimensions of the folded array
    folded_dims <- c(rep(n, p + 1), k)
    # Initialize the folded array
    folded <- array(0, dim = folded_dims)
    for (t in 1:k) {
        for (i in 1:n) {
            unfolded_slice <- unfolded[i, , t]
            folded_slice <- array(unfolded_slice, dim = rep(n, p))
            # Assign the reshaped slice to the folded array
            if (p == 1) {
                folded[i, , t] <- folded_slice
            } else if (p == 2) {
                folded[i, , , t] <- folded_slice
            } else if (p == 3) {
                folded[i, , , , t] <- folded_slice
            }
        }
    }
    return(folded)
}

unfold3d <- function(folded) {
    dims <- dim(folded)
    n <- dims[1]
    p <- length(dims) - 2
    k <- dim(folded)
    k <- k[length(k)]
    unfolded <- array(0, dim = c(n, n^p, k))
    for (t in 1:k) {
        for (i in 1:n) {
            # Extract the relevant slice from the folded array
            if (p == 1) {
                folded_slice <- folded[i, , t]
            } else if (p == 2) {
                folded_slice <- folded[i, , , t]
            } else if (p == 3) {
                folded_slice <- folded[i, , , , t]
            }
            # Flatten the folded slice and assign it to the unfolded array
            unfolded[i, , t] <- as.vector(folded_slice)
        }
    }
    return(unfolded)
}


fold4d <- function(unfolded, p) {
    dims <- dim(unfolded)
    nsim <- dims[length(dims)]
    h <- dims[length(dims) - 1]
    n <- dims[1]
    # Determine the dimensions of the folded array
    folded_dims <- c(rep(n, p + 1), h, nsim)
    # Initialize the folded array
    folded <- array(0, dim = folded_dims)
    for (sim in 1:nsim) {
        for (t in 1:h) {
            for (i in 1:n) {
                unfolded_slice <- unfolded[i, , t, sim]
                folded_slice <- array(unfolded_slice, dim = rep(n, p))
                # Assign the reshaped slice to the folded array
                if (p == 1) {
                    folded[i, , t, sim] <- folded_slice
                } else if (p == 2) {
                    folded[i, , , t, sim] <- folded_slice
                } else if (p == 3) {
                    folded[i, , , , t, sim] <- folded_slice
                }
            }
        }
    }
    return(folded)
}

unfold4d <- function(folded) {
    dims <- dim(folded)
    n <- dims[1]
    p <- length(dims) - 3
    dims <- dim(folded)
    k <- dims[length(dims) - 1]
    nsim <- dims[length(dims)]
    unfolded <- array(0, dim = c(n, n^p, k, nsim))
    for (sim in 1:nsim) {
        for (t in 1:k) {
            for (i in 1:n) {
                # Extract the relevant slice from the folded array
                if (p == 1) {
                    folded_slice <- folded[i, , t, sim]
                } else if (p == 2) {
                    folded_slice <- folded[i, , , t, sim]
                } else if (p == 3) {
                    folded_slice <- folded[i, , , , t, sim]
                }
                # Flatten the folded slice and assign it to the unfolded array
                unfolded[i, , t, sim] <- as.vector(folded_slice)
            }
        }
    }
    return(unfolded)
}


array4d_matrix_mean <- function(x)
{
    # assumes an s_series x n_series x h x nsim array
    n_series <- dim(x)[1]
    h <- dim(x)[3]
    nsim <- dim(x)[4]
    M <- array(0, dim = dim(x)[1:3])
    for (i in 1:h) {
        tmp <- x[,,i,]
        tmp <- apply(tmp, 3, identity, simplify = FALSE)
        M[,,i] <- Reduce(`+`, tmp)/nsim
    }
    return(M)
}

.gogarch_coskew_aggregate <- function(S, w)
{
    kronw <- kronecker(t(w), t(w))
    if (length(dim(S)) == 3) {
        n <- dim(S)[3]
        s <- rep(0, n)
        for (i in 1:n) {
            s[i] <- w %*% S[,,i] %*% kronw
        }
        return(s)
    } else {
        n <- dim(S)[3]
        sim <- dim(S)[4]
        s <- matrix(0, ncol = n, nrow = sim)
        for (i in 1:sim) {
            for (j in 1:n) {
                s[i,j] <- w %*% S[,,j,i] %*% kronw
            }
        }
        return(s)
    }
}

.gogarch_cokurt_aggregate <- function(K, w)
{
    kronw <- kronecker(t(w), kronecker(t(w), t(w)))
    if (length(dim(K)) == 3) {
        n <- dim(K)[3]
        k <- rep(0, n)
        for (i in 1:n) {
            k[i] <- w %*% K[,,i] %*% kronw
        }
        return(k)
    } else {
        n <- dim(K)[3]
        sim <- dim(K)[4]
        k <- matrix(0, ncol = n, nrow = sim)
        for (i in 1:sim) {
            for (j in 1:n) {
                k[i,j] <- w %*% K[,,j,i] %*% kronw
            }
        }
        return(k)
    }
}


ghyp_variance <- function(alpha, beta, delta, mu, lambda)
{
    zeta <- delta * sqrt(alpha^2 - beta^2)
    out <- besselK(zeta, lambda + 1)/(zeta * bessel_k(zeta, lambda)) + ((beta^2)/(alpha^2 - beta^2)) * ((bessel_k(zeta, lambda + 2)/bessel_k(zeta, lambda)) - ((bessel_k(zeta, lambda + 1)/bessel_k(zeta, lambda))^2))
    out <- delta^2 * out
    return(out)
}

ghyp_mu <- function(alpha, beta, delta, mu, lambda)
{
    zeta <- delta * sqrt(alpha^2 - beta^2)
    mu + (beta * delta)/(sqrt(alpha^2 - beta^2)) *  (besselK(zeta, lambda + 1)/(bessel_k(zeta, lambda)))
}

gh_support <- function(q, alpha, beta, delta, mu, lambda)
{
    zeta <- delta * sqrt(alpha^2 - beta^2)
    rho <- beta/alpha
    mu <- ghyp_mu(alpha, beta, delta, mu, lambda)
    sig <- sqrt(ghyp_variance(alpha, beta, delta, mu, lambda))
    return(qgh(q, mu = mu, sigma = sig, skew = rho, shape = zeta, lambda = lambda))

}
