#' GOGARCH Model specification
#'
#' @param y an xts matrix of pre-filtered (residuals) stationary data.
#' @param distribution a choice for the component distributions. Valid choices are
#' normal, normal inverse gaussian or generalized hyperbolic distribution.
#' @param model the GARCH model to use for each factor.
#' @param order the GARCH model order.
#' @param ica the Independent Component Analysis algorithm. Current only the
#' RADICAL algorithm is available.
#' @param components the number of components to extract in the pre-whitening
#' phase,
#' @param lambda_range for the generalized hyperbolic distribution, the range of
#' the lambda parameter.
#' @param shape_range for the generalized hyperbolic distribution, the range of
#' the shape parameter (zeta).
#' @param cond_mean an optional matrix of the conditional mean for the series.
#' @param ... additional arguments passed to the \code{\link{radical}} function.
#' @returns an object of class \dQuote{gogarch.spec}.
#' @export
gogarch_modelspec <- function(y, distribution = c("norm","nig","gh"), model = "garch", order = c(1,1), ica = "radical",
                              components = NCOL(y), lambda_range = c(-5, 5), shape_range = c(0.1, 25),
                              cond_mean = NULL, ...)
{

    distribution <- match.arg(distribution[1], c("norm","nig","gh"))
    ica <- match.arg(ica[1], "radical")
    spec <- list()
    if (!is.xts(y)) stop("\ny must be an xts matrix.")
    nobs <- NROW(y)
    n_series <- NCOL(y)
    series_names <- colnames(y)
    spec$target$y <- coredata(y)
    mu <- .cond_mean_spec(mu = cond_mean, n_series, NROW(y), series_names)
    spec$target$mu <- mu
    spec$target$index <- index(y)
    spec$target$sampling <-
    spec$garch$model <- model
    spec$garch$order <- order
    spec$garch$lambda_range <- lambda_range
    spec$garch$shape_range <- shape_range
    spec$ica$model <- ica
    spec$ica$components <- components
    spec$distribution <- distribution
    spec$series_names <- series_names
    spec$nobs <- nobs
    spec$n_series <- n_series
    class(spec) <- "gogarch.spec"
    return(spec)
}

.gogarch_estimate <- function(object, verbose = FALSE, ...)
{
    series <- estimate <- NULL
    # run ICA
    L <- list()
    ic <- radical(object$target$y, components = object$ica$components, demean = FALSE, trace = verbose, ...)
    components <- ic$S
    # run GARCH model on components
    garch_model <- .estimate_gogarch_components(components, object)
    # generate likelihood
    L$loglik <- .gogarch_log_likelihood(garch_model, ic$K)
    # extract covariance and correlation
    A <- ic$A
    V <- coredata(sigma(garch_model)^2)
    H <- .gogarch_covariance(V, A)
    R <- .gogarch_correlation(V, A)
    CS <- NULL
    CK <- NULL
    L$ica <- ic
    L$spec <- object
    L$univariate <- garch_model
    L$H <- H
    L$R <- R
    L$mu <- object$target$mu
    pmatrix <- lapply(garch_model, function(x) x$parmatrix)
    series_names <- paste0("ica_component.",1:length(garch_model))
    for (i in 1:length(pmatrix)) pmatrix[[i]][,series := series_names[i]]
    pmatrix <- rbindlist(pmatrix)
    L$parmatrix <- pmatrix
    # add back the number of parameters
    L$spec$npars <- NROW(pmatrix[estimate == 1]) + length(ic$U)
    class(L) <- "gogarch.estimate"
    return(L)
}

.garch_filter_model_ica <- function(object, y)
{
    m <- NCOL(y)
    new_fit <- NULL
    new_fit <- future_lapply(1:m, function(i) {
        tsfilter(object[[i]], y = y[,i])
    }, future.packages = "tsgarch", future.seed = TRUE)
    new_fit <- eval(new_fit)
    names(new_fit) <- names(object)
    new_fit <- to_multi_estimate(new_fit)
    return(new_fit)
}

.gogarch_filter <- function(object, y = NULL, cond_mean = NULL, ...)
{
    elapsed <- Sys.time()
    if (!is.xts(y)) stop("\ny must be an xts object.")
    if (!is.null(y)) {
        mu <- .cond_mean_spec(mu = cond_mean, object$spec$n_series, NROW(y), object$spec$series_names)
        object$spec$target$mu <- rbind(object$spec$target$mu, mu)
        object$mu <- rbind(object$mu, mu)
    }
    y <- .check_y_filter(object, y = y)
    new_S <- filter_radical(y, demean = FALSE, object$ica$W, mu = NULL)
    S <- rbind(object$ica$S, new_S)
    # filter y
    if (is.null(object$spec$target$original_size)) {
        object$spec$target$original_size <- NROW(object$spec$target$y)
    }
    garch_model <- .garch_filter_model_ica(object$univariate, xts(new_S, index(y)))
    object$loglik <- .gogarch_log_likelihood(garch_model, object$ic$K)
    object$spec$nobs <- NROW(S)
    object$spec$target$y <- rbind(object$spec$target$y, coredata(y))
    object$spec$target$index <- c(object$spec$target$index, index(y))
    V <- coredata(sigma(garch_model)^2)
    H <- .gogarch_covariance(V, object$ica$A)
    R <- .gogarch_correlation(V, object$ica$A)
    CS <- NULL
    CK <- NULL
    object$ica$S <- S
    object$univariate <- garch_model
    object$H <- H
    object$R <- R
    class(object) <- "gogarch.estimate"
    return(object)
}


.gogarch_predict <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"),
                             forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
{
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    if (is.null(forc_dates)) {
        forc_dates <- .forecast_dates(forc_dates, h = h, sampling = object$univariate[[1]]$spec$target$sampling,
                                      last_index = tail(object$spec$target$index,1))
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h.")
    }
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    sim_method <- match.arg(sim_method, c("parametric", "bootstrap"))
    p <- lapply(1:length(object$univariate), function(i) {
        pred <- predict(object$univariate[[i]], h = h, nsim = nsim, sim_method = sim_method, forc_dates = forc_dates, seed = seed)
        # tsgarch returns NULL when h = 1 for sigma_sim (there is no uncertainty). But to be consistent with the approach for dcc
        # cgarch, we return an array
        if (is.null(pred$sigma_sim) & h == 1) {
            sim_sigma <- matrix(as.numeric(pred$sigma), ncol = h, nrow = nsim)
            colnames(sim_sigma) <- as.character(forc_dates)
            class(sim_sigma) <- "tsmodel.distribution"
            pred$sigma_sim <- sim_sigma
            return(pred)
        } else {
            return(pred)
        }
    })
    # need conditional mean
    sim_mu <- lapply(p, function(x) x$distribution)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, length(p)))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    # translate back
    sim_mu <- apply(sim_mu, 3, FUN = function(x) x %*% object$ica$A, simplify = FALSE)
    sim_mu <- array(unlist(sim_mu), dim = c(h, object$spec$n_series, nsim))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
    elapsed <- Sys.time() - elapsed
    out <- list(univariate = p, mu = sim_mu, cond_mean = mu, ica = object$ic, parmatrix = object$parmatrix, spec = object$spec, forc_dates = forc_dates,
                h = h, nsim = nsim, n_series = object$spec$n_series, series_names = object$spec$series_names, elapsed = elapsed)
    class(out) <- "gogarch.predict"
    return(out)
}

.gogarch_simulate <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0, cond_mean = NULL, sim_method = c("parametric", "bootstrap"), ...)
{
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    sim_method <- match.arg(sim_method[1], c("parametric", "bootstrap"))
    sim <- lapply(1:length(object$univariate), function(i){
        xspec <- object$univariate[[i]]$spec
        xspec$parmatrix <- copy(object$univariate[[i]]$parmatrix)
        simulate(xspec, h = h, nsim = nsim, sim_method = sim_method, burn = burn, seed = seed)
    })
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    sim_mu <- lapply(sim, function(x) x$series)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, length(sim)))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    # translate back
    sim_mu <- apply(sim_mu, 3, FUN = function(x) x %*% object$ica$A, simplify = FALSE)
    sim_mu <- array(unlist(sim_mu), dim = c(h, object$spec$n_series, nsim))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
    elapsed <- Sys.time() - elapsed
    out <- list(univariate = sim, mu = sim_mu, ica = object$ic, parmatrix = object$parmatrix, spec = object$spec,
                h = h, nsim = nsim, elapsed = elapsed)
    class(out) <- "gogarch.simulate"
    return(out)
}

# always return [h n_series nsim]
# NOTE: need to deal with dimensionality reduced system
.gogarch_fitted <- function(object)
{
    if (inherits(object, "gogarch.estimate")) {
        # strictly speaking this is not requires since the fitted from the IC is always zero mean
        x <- xts(object$mu, object$spec$target$index)
        colnames(x) <- object$spec$series_names
    } else if (inherits(object, "gogarch.predict")) {
        x <- object$mu
        attr(x,"index") <- as.character(object$spec$target$index)
        attr(x,"series") <- object$spec$series_names
    } else if (inherits(object, "gogarch.simulate")) {
        x <- object$mu
        attr(x,"series") <- object$spec$series_names
    }
    return(x)
}

.gogarch_dskewness <- function(object)
{
    parameter <- NULL
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    sk <- dskewness(object$spec$distribution, skew = skew, shape = shape, lambda = lambda)
    return(sk)
}


.gogarch_dkurtosis <- function(object)
{
    parameter <- NULL
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    ku <- dkurtosis(object$spec$distribution, skew = skew, shape = shape, lambda = lambda) + 3
    return(ku)
}


.gogarch_port_skewness_estimate <- function(object, weights, sigma)
{
    A <- object$ica$A
    sig <- coredata(sigma(object$univariate))
    if (object$spec$distribution != 'norm') {
        sk <- .gogarch_dskewness(object)
        sk <- matrix(sk, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^3
        S <- .gogarch_skewness_weighted(A, sk, weights)/sigma^3
        return(S)
    } else {
        warning("\nmultivariate normal does not have skewness")
    }
}

.gogarch_port_kurtosis_estimate <- function(object, weights, sigma)
{
    A <- object$ica$A
    sig <- coredata(sigma(object$univariate))
    if (object$spec$distribution != 'norm') {
        ku <- .gogarch_dkurtosis(object)
        ku <- matrix(ku, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^4
        K <- .gogarch_kurtosis_weighted(A, K = ku, V = sig^2, w = weights)/sigma^4
        return(K)
    } else {
        warning("\nmultivariate normal does not have excess kurtosis")
    }
}


.gogarch_port_sigma_simulate <- function(object, weights)
{
    A <- object$ica$A
    m <- NROW(A)
    n <- object$h
    nsim <- object$nsim
    if (is(object, "gogarch.predict")) {
        sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma_sim)
    } else if (is(object, "gogarch.simulate")) {
        sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma)
    }
    sig <- array(unlist(sig), dim = c(nsim, n, m))
    sig <- aperm(sig, perm = c(2,3,1))
    S <- matrix(0, nrow = nsim, ncol = n)
    for (i in 1:nsim) {
        sig_mat <- .retain_dimensions_array(sig, i)
        S[i,] <- .gogarch_covariance_weighted(sig_mat^2, A, weights)
    }
    S <- sqrt(S)
    return(S)
}

.gogarch_port_skewness_simulate <- function(object, weights, sigma)
{
    A <- object$ica$A
    m <- NROW(A)
    n <- object$h
    nsim <- object$nsim
    if (is(object, "gogarch.predict")) {
        sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma_sim)
    } else if (is(object, "gogarch.simulate")) {
        sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma)
    }
    sig <- array(unlist(sig), dim = c(nsim, n, m))
    sig <- aperm(sig, perm = c(2,3,1))
    if (object$spec$distribution != 'norm') {
        sk <- .gogarch_dskewness(object)
        S <- matrix(0, nrow = nsim, ncol = n)
        for (i in 1:nsim) {
            sig_mat <- .retain_dimensions_array(sig, i)
            sksim <- matrix(sk, ncol = m, nrow = n, byrow = TRUE) * sig_mat^3
            S[i,] <- .gogarch_skewness_weighted(A, sksim, weights)
        }
        S <- S/sigma^3
        return(S)
    } else {
        warning("\nmultivariate normal does not have skewness")
    }
}

.gogarch_port_kurtosis_simulate <- function(object, weights, sigma)
{
    A <- object$ica$A
    m <- NROW(A)
    n <- object$h
    nsim <- object$nsim
    if (is(object, "gogarch.predict")) {
        sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma_sim)
    } else if (is(object, "gogarch.simulate")) {
        sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma)
    }
    sig <- array(unlist(sig), dim = c(nsim, n, m))
    sig <- aperm(sig, perm = c(2,3,1))
    if (object$spec$distribution != 'norm') {
        ku <- .gogarch_dkurtosis(object)
        kusim <- matrix(ku, ncol = m, nrow = n, byrow = TRUE)
        K <- .gogarch_kurtosis_weighted_sim(A, sig, kusim, weights, nsim, n)
        K <- K/sigma^4
        return(K)
    } else {
        warning("\nmultivariate normal does not have kurtosis")
    }
}


# convolution of distributions given weights
# d, p, q, r of convoluted object



.nig_fft <- function(w_pars, fft_support, fft_step, fft_by)
{
    n <- dim(w_pars)[3]
    out <- future_lapply(1:n, function(i) {
        if (is.null(fft_support)) {
            x_min <- min(apply(cbind(w_pars[,1,i], w_pars[,2,i], w_pars[,3,i], w_pars[,4,i]), 1,
                              FUN = function(x) .qnig(0.0000001, alpha = x[1], beta = x[2], delta = x[3], mu = x[4])))
            x_max <- max(apply(cbind(w_pars[,1,i], w_pars[,2,i], w_pars[,3,i], w_pars[,4,i]), 1,
                              FUN = function(x) .qnig(1 - 0.0000001, alpha = x[1], beta = x[2], delta = x[3], mu = x[4])))
            zz <- seq(x_min, x_max, by = fft_by)
        } else {
            x_min <- fft_support[1]
            x_max <- fft_support[2]
            zz <- seq(x_min, x_max, by = fft_by)
        }
        pdf <- .cfinvnig(z = zz, step = fft_step, alpha = w_pars[,1,i], beta = w_pars[,2,i], delta = w_pars[,3,i], mu = w_pars[,4,i])
        attr(pdf, "support") <- c(x_min, x_max)
        return(pdf)
    }, future.seed = TRUE)
    return(out)
}

.gh_fft <- function(w_pars, fft_support, fft_step, fft_by)
{
    n <- dim(w_pars)[3]
    out <- future_lapply(1:n, function(i) {
        if (is.null(fft_support)) {
            tmp <- cbind(w_pars[,1,i], w_pars[,2,i], w_pars[,3,i], w_pars[,4,i], w_pars[,5,i])
            x_min <- min(sapply(1:NROW(tmp), function(j) {
                ans <- try(gh_support(0.00000001, alpha = tmp[j,1], beta = tmp[j,2], delta = tmp[j,3], mu = tmp[j,4], lambda = tmp[j,5]), silent = TRUE)
                if (inherits(ans, 'try-error')) return(-1) else return(ans)
            }), na.rm = TRUE)
            x_max <- max(sapply(1:NROW(tmp), function(j) {
                ans <- try(gh_support(1 - 0.00000001, alpha = tmp[j,1], beta = tmp[j,2], delta = tmp[j,3], mu = tmp[j,4], lambda = tmp[j,5]), silent = TRUE)
                if (inherits(ans, 'try-error')) return(1) else return(ans)
                }), na.rm = TRUE)
            zz <- seq(x_min, x_max, by = fft_by)
        } else {
            x_min <- fft_support[1]
            x_max <- fft_support[2]
            zz <- seq(x_min, x_max, by = fft_by)
        }
        pdf <- .cfinvghyp(z = zz, step = fft_step, alpha = w_pars[,1,i], beta = w_pars[,2,i], delta = w_pars[,3,i], mu = w_pars[,4,i], lambda = w_pars[,5,i])
        attr(pdf, "support") <- c(x_min, x_max)
        return(pdf)
    }, future.seed = TRUE)
    return(out)
}


.gh_fft_sim <- function(w_pars, fft_support, fft_step, fft_by)
{
    n <- dim(w_pars)[3]
    out <- lapply(1:n, function(i) {
        if (is.null(fft_support)) {
            tmp <- cbind(w_pars[,1,i], w_pars[,2,i], w_pars[,3,i], w_pars[,4,i], w_pars[,5,i])
            x_min <- min(sapply(1:NROW(tmp), function(j) {
                ans <- try(gh_support(0.00000001, alpha = tmp[j,1], beta = tmp[j,2], delta = tmp[j,3], mu = tmp[j,4], lambda = tmp[j,5]), silent = TRUE)
                if (inherits(ans, 'try-error')) return(-1) else return(ans)
            }), na.rm = TRUE)
            x_max <- max(sapply(1:NROW(tmp), function(j) {
                ans <- try(gh_support(1 - 0.00000001, alpha = tmp[j,1], beta = tmp[j,2], delta = tmp[j,3], mu = tmp[j,4], lambda = tmp[j,5]), silent = TRUE)
                if (inherits(ans, 'try-error')) return(1) else return(ans)
            }), na.rm = TRUE)
            zz <- seq(x_min, x_max, by = fft_by)
        } else {
            x_min <- fft_support[1]
            x_max <- fft_support[2]
            zz <- seq(x_min, x_max, by = fft_by)
        }
        pdf <- .cfinvghyp(z = zz, step = fft_step, alpha = w_pars[,1,i], beta = w_pars[,2,i], delta = w_pars[,3,i], mu = w_pars[,4,i], lambda = w_pars[,5,i])
        attr(pdf, "support") <- c(x_min, x_max)
        return(pdf)
    })
    return(out)
}

.gogarch_newsimpact_covariance <- function(object, pair = c(1,2), factor = c(1,2), ...)
{
    # factors
    A <- object$ica$A
    m <- NROW(A)
    newsimpact_list <- vector(mode = "list", length = m)
    pair_names <- object$spec$series_names[pair]
    newsimpact_list[[1]] <- newsimpact(object = object$univariate[[1]], epsilon = NULL)
    # we want a common epsilon
    z_eps <- newsimpact_list[[1]]$x
    for (i in 2:m) newsimpact_list[[i]] <- newsimpact(epsilon = z_eps, object = object$univariate[[i]])
    sigmas <- sapply(newsimpact_list, FUN = function(x) x$y)
    N <- NROW(sigmas)
    H <- array(NA, dim = c(NCOL(A), NCOL(A), N))
    ni <- matrix(NA, ncol = N, nrow = N)
    # find the no-shock index (z_eps = 0)
    zero_shock <- as.numeric(which(round(newsimpact_list[[factor[1]]]$x,12) == 0)[1])
    # to get the off diagonals we need to keep 1 fixed and
    # revolve the epsilon of the other
    s2 <- sapply(newsimpact_list, FUN = function(x) x$y)
    rownames(s2) <- round(as.numeric(z_eps),4)
    np <- 1:m
    for (j in 1:N) {
        for (i in 1:N) {
            sigmas <- s2[zero_shock,]
            sigmas1 <- as.numeric(s2[j, factor[1]])
            sigmas2 <- as.numeric(s2[i, factor[2]])
            sigmas[factor[1]] <- sigmas1
            sigmas[factor[2]] <- sigmas2
            H <- t(A) %*% diag(sigmas) %*% A
            ni[j,i] <- H[pair[1], pair[2]]
        }
    }
    out <- list(nisurface = ni, axis = z_eps, pairs = pair, pair_names = pair_names, factors = factor, type = "covariance", model = "gogarch", dynamics = object$spec$garch$model, A = A)
    class(out) <- c("tsmarch.newsimpact")
    return(out)
}


.gogarch_newsimpact_correlation <- function(object, pair = c(1,2), factor = c(1,2), ...)
{
    # factors
    A <- object$ica$A
    m <- NROW(A)
    newsimpact_list <- vector(mode = "list", length = m)
    pair_names <- object$spec$series_names[pair]
    newsimpact_list[[1]] <- newsimpact(object = object$univariate[[1]], epsilon = NULL)
    # we want a common epsilon
    z_eps <- newsimpact_list[[1]]$x
    for (i in 2:m) newsimpact_list[[i]] <- newsimpact(epsilon = z_eps, object = object$univariate[[i]])
    sigmas <- sapply(newsimpact_list, FUN = function(x) x$y)
    N <- NROW(sigmas)
    H <- array(NA, dim = c(NCOL(A), NCOL(A), N))
    ni <- matrix(NA, ncol = N, nrow = N)
    # find the no-shock index (z_eps = 0)
    zero_shock <- as.numeric(which(round(newsimpact_list[[factor[1]]]$x,12) == 0)[1])
    # to get the off diagonals we need to keep 1 fixed and
    # revolve the epsilon of the other
    s2 <- sapply(newsimpact_list, FUN = function(x) x$y)
    rownames(s2) <- round(as.numeric(z_eps),4)
    np <- 1:m
    for (j in 1:N) {
        for (i in 1:N) {
            sigmas <- s2[zero_shock,]
            sigmas1 <- as.numeric(s2[j, factor[1]])
            sigmas2 <- as.numeric(s2[i, factor[2]])
            sigmas[factor[1]] <- sigmas1
            sigmas[factor[2]] <- sigmas2
            H <- t(A) %*% diag(sigmas) %*% A
            R <- cov2cor(H)
            ni[j,i] <- R[pair[1], pair[2]]
        }
    }
    out <- list(nisurface = ni, axis = z_eps, pairs = pair, pair_names = pair_names, factors = factor, type = "correlation", model = "gogarch", dynamics = object$spec$garch$model, A = A)
    class(out) <- c("tsmarch.newsimpact")
    return(out)
}

.gogarch_news_impact_surface <- function(x, ...)
{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    heading <- paste0(toupper(x$model)," [",toupper(x$dynamics),"] News Impact ", upper_case_i(x$type,1,1)," Surface\n",x$pair_names[1],"[",x$factors[1],"]:",x$pair_names[2],"[",x$factors[2],"]")
    layout(matrix(c(1,1,2,1,1,3), 2, 3, byrow = TRUE))
    factors <- x$factors
    axis_names <- c(paste0("F_",factors[1]," [t-1]"), paste0("F_",factors[2]," [t-1]"))
    x1 <- drapecol(x$nisurface, col = topo.colors(100), NAcol = "white")
    par(mar = c(1.5,1.5,2,1.5))
    persp(  x = x$axis,
            y = x$axis,
            z = x$nisurface,  col = x1, theta = 45, phi = 25, expand = 0.5,
            ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = axis_names[1],
            ylab = axis_names[2], zlab = "",
            cex.axis = 0.7, cex.main = 0.8, main = heading)
    par(mar = c(2,2,2,2))
    cols_factor1 <- c("cadetblue","steelblue")
    barplot(x$A[factors,x$pairs[1]], col = cols_factor1, main = paste0(x$pair_names[1]," Factor Loadings"), cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, names.arg = paste0("F_",factors),  cex.main = 0.7)
    abline(h = 0)
    box()
    cols_factor2 <- c("cadetblue","steelblue")
    barplot(x$A[factors,x$pairs[2]], col = cols_factor2, main = paste0(x$pair_names[2]," Factor Loadings"), cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, cex.main = 0.7)
    abline(h = 0)
    box()
    return(invisible(x))
}

.gogarch_convolution_prediction_d <- function(object,  weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1))
{
    elapsed <- Sys.time()
    n <- object$h
    nsim <- object$nsim
    A <- object$ica$A
    m <- NROW(A)
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    parmatrix <- .get_garch_parameters(object)
    # transform to alpha, beta, delta, mu, lambda
    if (object$spec$distribution == "nig") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            nigtransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i])
        }))
        # sort to have alpha beta delta mu
        pmatrix_std <- pmatrix_std[,c(4,3,2,1)]
    } else if (object$spec$distribution == "gh") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            ghyptransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i], lambda = parmatrix[3,i])
        }))
        pmatrix_std <- cbind(pmatrix_std, parmatrix[3,])
        pmatrix_std <- pmatrix_std[,c(4,3,2,1,5)]
        colnames(pmatrix_std)[5] <- "lambda"
    } else {
        stop("\nconvolution not required for normal distribution.")
    }
    # no uncertainty for h = 1
    sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma_sim)
    sig <- array(unlist(sig), dim = c(nsim, n, m))
    sig <- aperm(sig, perm = c(2,3,1))
    mu <- object$mu
    w_mu <- matrix(NA, ncol = n, nrow = nsim)
    for (i in 1:n) {
        w_mu[,i] <- t(mu[i,,]) %*% w[i,]
    }
    Q <- future_lapply(1:nsim, function(s) {
        w_hat <- matrix(NA, ncol = m, nrow = n)
        for (i in 1:n) {
            dS <- diag(sig[i,,s])
            w_hat[i,] <- w[i,] %*% (t(A) %*% dS)
        }
        if (object$spec$distribution == "nig") {
            w_pars <- array(data = NA, dim = c(m, 4, n))
            # Scaling of parameters (Blaesild)
            for (j in 1:m) {
                tmp <- matrix(pmatrix_std[j,1:4], ncol = 4, nrow = n, byrow = TRUE)
                tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j])
                w_pars[j,,] <- as.array(t(tmp), dim = c(1, 4, n))
            }
            out <- .nig_fft(w_pars, fft_support, fft_step, fft_by)
        } else if (object$spec$distribution == "gh") {
            w_pars <- array(data = NA, dim = c(m, 5, n))
            # Scaling of parameters (Blaesild)
            for (j in 1:m) {
                tmp <- matrix(pmatrix_std[j,1:5], ncol = 5, nrow = n, byrow = TRUE)
                tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j], rep(1, n))
                w_pars[j,,] <- as.array(t(tmp), dim = c(1, 5, n))
            }
            out <- .gh_fft_sim(w_pars, fft_support, fft_step, fft_by)
        }
    }, future.seed = TRUE)
    elapsed <- Sys.time() - elapsed
    L <- list(y = Q, distribution = object$spec$distribution, mu = w_mu, fft_step = fft_step, fft_by = fft_by, elapsed = elapsed)
    class(L) <- "gogarch.fftsim"
    return(L)
}


.gogarch_convolution_prediction <- function(object,  weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1))
{
    elapsed <- Sys.time()
    n <- object$h
    A <- object$ica$A
    m <- NROW(A)
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    parmatrix <- .get_garch_parameters(object)
    # transform to alpha, beta, delta, mu, lambda
    if (object$spec$distribution == "nig") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            nigtransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i])
        }))
        # sort to have alpha beta delta mu
        pmatrix_std <- pmatrix_std[,c(4,3,2,1)]
    } else if (object$spec$distribution == "gh") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            ghyptransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i], lambda = parmatrix[3,i])
        }))
        pmatrix_std <- cbind(pmatrix_std, parmatrix[3,])
        pmatrix_std <- pmatrix_std[,c(4,3,2,1,5)]
        colnames(pmatrix_std)[5] <- "lambda"
    } else {
        stop("\nconvolution not required for normal distribution.")
    }
    # no uncertainty for h = 1
    sig <- do.call(cbind, lapply(object$univariate, function(x) coredata(x$sigma)))
    # this is the same as object$cond_mean when present since the distribution has been
    # recentered around this values.
    mu <- t(apply(object$mu, 1, rowMeans))
    w_hat <- matrix(NA, ncol = m, nrow = n)
    w_mu <- rep(NA, n)
    for (i in 1:n) {
        dS <- diag(sig[i, ])
        w_hat[i,] <- w[i,] %*% (t(A) %*% dS)
        w_mu[i] <- mu[i, ] %*% w[i, ]
    }
    if (object$spec$distribution == "nig") {
        w_pars <- array(data = NA, dim = c(m, 4, n))
        # Scaling of parameters (Blaesild)
        for (j in 1:m) {
            tmp <- matrix(pmatrix_std[j,1:4], ncol = 4, nrow = n, byrow = TRUE)
            tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j])
            w_pars[j,,] <- as.array(t(tmp), dim = c(1, 4, n))
        }
        out <- .nig_fft(w_pars, fft_support, fft_step, fft_by)
    } else if (object$spec$distribution == "gh") {
        w_pars <- array(data = NA, dim = c(m, 5, n))
        # Scaling of parameters (Blaesild)
        for (j in 1:m) {
            tmp <- matrix(pmatrix_std[j,1:5], ncol = 5, nrow = n, byrow = TRUE)
            tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j], rep(1, n))
            w_pars[j,,] <- as.array(t(tmp), dim = c(1, 5, n))
        }
        out <- .gh_fft(w_pars, fft_support, fft_step, fft_by)
    }
    L <- list(y = out, distribution = object$spec$distribution, mu = w_mu, fft_step = fft_step, fft_by = fft_by)
    class(L) <- "gogarch.fft"
    return(L)
}

.gogarch_convolution_simulation_d <- function(object,  weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1))
{
    elapsed <- Sys.time()
    n <- object$h
    nsim <- object$nsim
    A <- object$ica$A
    m <- NROW(A)
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    parmatrix <- .get_garch_parameters(object)
    # transform to alpha, beta, delta, mu, lambda
    if (object$spec$distribution == "nig") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            nigtransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i])
        }))
        # sort to have alpha beta delta mu
        pmatrix_std <- pmatrix_std[,c(4,3,2,1)]
    } else if (object$spec$distribution == "gh") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            ghyptransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i], lambda = parmatrix[3,i])
        }))
        pmatrix_std <- cbind(pmatrix_std, parmatrix[3,])
        pmatrix_std <- pmatrix_std[,c(4,3,2,1,5)]
        colnames(pmatrix_std)[5] <- "lambda"
    } else {
        stop("\nconvolution not required for normal distribution.")
    }
    # no uncertainty for h = 1
    sig <- lapply(1:m, function(i) object$univariate[[i]]$sigma)
    sig <- array(unlist(sig), dim = c(nsim, n, m))
    sig <- aperm(sig, perm = c(2,3,1))
    mu <- object$mu
    w_mu <- matrix(NA, ncol = n, nrow = nsim)
    for (i in 1:n) {
        w_mu[,i] <- t(mu[i,,]) %*% w[i,]
    }
    Q <- future_lapply(1:nsim, function(s) {
        w_hat <- matrix(NA, ncol = m, nrow = n)
        for (i in 1:n) {
            dS <- diag(sig[i,,s])
            w_hat[i,] <- w[i,] %*% (t(A) %*% dS)
        }
        if (object$spec$distribution == "nig") {
            w_pars <- array(data = NA, dim = c(m, 4, n))
            # Scaling of parameters (Blaesild)
            for (j in 1:m) {
                tmp <- matrix(pmatrix_std[j,1:4], ncol = 4, nrow = n, byrow = TRUE)
                tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j])
                w_pars[j,,] <- as.array(t(tmp), dim = c(1, 4, n))
            }
            out <- .nig_fft(w_pars, fft_support, fft_step, fft_by)
        } else if (object$spec$distribution == "gh") {
            w_pars <- array(data = NA, dim = c(m, 5, n))
            # Scaling of parameters (Blaesild)
            for (j in 1:m) {
                tmp <- matrix(pmatrix_std[j,1:5], ncol = 5, nrow = n, byrow = TRUE)
                tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j], rep(1, n))
                w_pars[j,,] <- as.array(t(tmp), dim = c(1, 5, n))
            }
            out <- .gh_fft_sim(w_pars, fft_support, fft_step, fft_by)
        }
    }, future.seed = TRUE)
    elapsed <- Sys.time() - elapsed
    L <- list(y = Q, distribution = object$spec$distribution, mu = w_mu, fft_step = fft_step, fft_by = fft_by, elapsed = elapsed)
    class(L) <- "gogarch.fftsim"
    return(L)
}


.gogarch_convolution_simulation <- function(object,  weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1))
{
    elapsed <- Sys.time()
    n <- object$h
    A <- object$ica$A
    m <- NROW(A)
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    parmatrix <- .get_garch_parameters(object)
    # transform to alpha, beta, delta, mu, lambda
    if (object$spec$distribution == "nig") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            nigtransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i])
        }))
        # sort to have alpha beta delta mu
        pmatrix_std <- pmatrix_std[,c(4,3,2,1)]
    } else if (object$spec$distribution == "gh") {
        pmatrix_std <- do.call(rbind, lapply(1:ncol(parmatrix), function(i){
            ghyptransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i], lambda = parmatrix[3,i])
        }))
        pmatrix_std <- cbind(pmatrix_std, parmatrix[3,])
        pmatrix_std <- pmatrix_std[,c(4,3,2,1,5)]
        colnames(pmatrix_std)[5] <- "lambda"
    } else {
        stop("\nconvolution not required for normal distribution.")
    }
    # no uncertainty for h = 1
    sig <- do.call(cbind, lapply(1:m, function(i) sqrt(colMeans(object$univariate[[i]]$sigma^2))))
    mu <- t(apply(object$mu, 1, rowMeans))
    w_hat <- matrix(NA, ncol = m, nrow = n)
    w_mu <- rep(NA, n)
    for (i in 1:n) {
        dS <- diag(sig[i, ])
        w_hat[i,] <- w[i,] %*% (t(A) %*% dS)
        w_mu[i] <- mu[i, ] %*% w[i, ]
    }
    if (object$spec$distribution == "nig") {
        w_pars <- array(data = NA, dim = c(m, 4, n))
        # Scaling of parameters (Blaesild)
        for (j in 1:m) {
            tmp <- matrix(pmatrix_std[j,1:4], ncol = 4, nrow = n, byrow = TRUE)
            tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j])
            w_pars[j,,] <- as.array(t(tmp), dim = c(1, 4, n))
        }
        out <- .nig_fft(w_pars, fft_support, fft_step, fft_by)
    } else if (object$spec$distribution == "gh") {
        w_pars <- array(data = NA, dim = c(m, 5, n))
        # Scaling of parameters (Blaesild)
        for (j in 1:m) {
            tmp <- matrix(pmatrix_std[j,1:5], ncol = 5, nrow = n, byrow = TRUE)
            tmp <- tmp * cbind(1/abs(w_hat[1:n, j]), 1/w_hat[1:n,j], abs(w_hat[1:n,j]), w_hat[1:n,j], rep(1, n))
            w_pars[j,,] <- as.array(t(tmp), dim = c(1, 5, n))
        }
        out <- .gh_fft(w_pars, fft_support, fft_step, fft_by)
    }
    L <- list(y = out, distribution = object$spec$distribution, mu = w_mu, fft_step = fft_step, fft_by = fft_by)
    class(L) <- "gogarch.fft"
    return(L)
}

