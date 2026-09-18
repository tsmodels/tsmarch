.calendar_eom <- function(date, ...)
{
    if (!is(date, "Date")) date <- as.Date(date)
    # Add a month, then subtract a day:
    date.lt <- as.POSIXlt(date, format = "%Y-%m-%d", tz = tz(date))
    mon <- date.lt$mon + 2
    year <- date.lt$year
    # If month was December add a year
    year <- year + as.integer(mon == 13)
    mon[mon == 13] <- 1
    iso <- ISOdate(1900 + year, mon, 1, hour = 0, tz = tz(date))
    result <- as.POSIXct(iso) - 86400 # subtract one day
    result <- result + (as.POSIXlt(iso)$isdst - as.POSIXlt(result)$isdst)*3600
    result <- as.Date(result)
    return(result)
}

.future_dates <- function(start, frequency, n = 1)
{
    if (frequency %in% c("days", "weeks", "months","years")) {
        switch(frequency,
               "days"   = as.Date(start) %m+% days(1:n),
               "weeks"  = as.Date(start) %m+% weeks(1:n),
               "months" = .calendar_eom(as.Date(start) %m+% months(1:n)),
               "years"  = as.Date(start) %m+% years(1:n))
    } else if (grepl("secs|mins|hours|",frequency)) {
        # Add one extra point and eliminate first one
        seq(as.POSIXct(start), length.out = n + 1, by = frequency)[-1]
    } else{
        as.Date(start) + (1:n)
    }
}

.forecast_dates <- function(forc_dates = NULL, h = 1, sampling, last_index)
{
    if (is.null(forc_dates)) {
        forc_dates = .future_dates(last_index, frequency = sampling, n = h)
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h")
        if (any(forc_dates <= last_index)) stop("\nforc_dates must be stricly greater than in-sample dates and times.")
    }
    return(forc_dates)
}

.sign_matrix <- function(x)
{
    ans = (-sign(x) + 1)/2
    ans[ans == 0.5] = 0
    ans
}

.repmat <- function(a, n, m)
{
    kronecker(matrix(1, n, m), a)
}


.lower_tri_matrix <- function(x, diag = FALSE)
{
    dims <- dim(x)
    if (length(dims) == 2) {
        out <- x[lower.tri(x, diag = diag)]
    } else if (length(dims) == 3) {
        out <- lapply(1:dims[3], function(i) x[,,i][lower.tri(x[,,i], diag = diag)])
        out <- do.call(rbind, out)
    } else if (length(dims) == 4) {
        n <- .tri_dimension(dims[1], diag)
        out <- array(0, dim = c(dims[3], n, dims[4]))
        for (i in 1:dims[4]) {
            tmp <- lapply(1:dims[3], function(j) x[,,j,i][lower.tri(x[,,j,i], diag = diag)])
            tmp <- do.call(rbind, tmp)
            out[,,i] <- tmp
        }
    }
    return(out)
}

.compose_tri_matrix <- function(x, dims, diag = FALSE) {
    if (length(dims) == 2) {
        if (!diag) {
            out <- matrix(0, dims[1], dims[2])
            out[lower.tri(out, diag = FALSE)] <- x
            out <- out + t(out)
            diag(out) <- 1
        } else {
            out <- matrix(0, dims[1], dims[2])
            out[lower.tri(out, diag = FALSE)] <- x
            out <- out + t(out)
            diag(out) <- diag(out)/2
        }
    } else if (length(dims) == 3) {
        out <- array(0, dims)
        if (!diag) {
            for (i in 1:dims[3]) {
                tmp <- matrix(0, dims[1], dims[2])
                tmp[lower.tri(tmp, diag = FALSE)] <- x[i,]
                tmp <- tmp + t(tmp)
                diag(tmp) <- 1
                out[,,i] <- tmp
            }
        } else {
            for (i in 1:dims[3]) {
                tmp <- matrix(0, dims[1], dims[2])
                tmp[lower.tri(tmp, diag = TRUE)] <- x[i,]
                tmp <- tmp + t(tmp)
                diag(tmp) <- diag(tmp)/2
                out[,,i] <- tmp
            }
        }
    } else {
        stop("\nonly 2d and 3d objects currently supported.")
    }
    return(out)
}


.lower_tri_dimension <- function(n, diag = TRUE)
{
    if (diag) {
        # positive solution to the inverse of (x^2 + x)/2
        return(-1/2 + sqrt(1 + 8*n)/2)
    } else {
        # positive solution to the inverse of (x^2 - x)/2
        return(1/2 + sqrt(1 + 8*n)/2)
    }
}

.tri_dimension <- function(n, diag = TRUE)
{
    if (diag) {
        return((n^2 + n)/2)
    } else {
        return((n^2 - n)/2)
    }
}

.lower_tri_pairs_index <- function(pair, n, diag = TRUE)
{
    # n is the vector size
    # convert to matrix dimension
    # solution by alexiosg
    n <- .lower_tri_dimension(n, diag = diag)
    if (max(pair) > n) stop("\nmax(pair)>matrix dimensions.")
    if (min(pair) <= 0) stop("\nmin(pair) <= 0.")
    # use decreasing order for lower triangular matrix
    pair <- sort(pair, decreasing = TRUE)
    tmp <- matrix(0, n, n)
    idx <- which(lower.tri(tmp, diag = diag), arr.ind = TRUE)
    # return pair combination index
    which(idx[,1] == pair[1] & idx[,2] == pair[2])
}

.is_equal <- function(p0, p1, tol = min(1e-05/2, .Machine$double.eps^.7)) {
    return(abs(p0 - p1) < tol)
}

.is_inside <- function(p0, pmat, tol = min(1e-05/2,.Machine$double.eps^.7))
{
    list1 <- lapply(1:nrow(pmat), function(x) {
        (p0 + tol > pmat[x, 1]) & (p0 - tol < pmat[x,2])
        })
    out <- apply(matrix(unlist(list1), ncol = nrow(pmat)), 1, any)
    return(out)
}

.is_equal_01 <- function(x) {
    return(.is_equal(x, 0) | .is_equal(x, 1))
}

.set_equal <- function(x, y, tol = 1e-7) {
    x1 <- round(2 * x / tol, 0)
    y1 <- round(2 * y / tol, 0)
    z <- x
    m <- match(x1, y1)
    valid_m <- !is.na(m)
    z[valid_m] <- y[m[valid_m]]
    return(z)
}

upper_case_i <- function(x, start = 1, end = 1) {
    substr(x, start, end) <- toupper(substr(x, start, end))
    x
}

bdiag <- function(...)
{
    if (nargs() == 1) x <- as.list(...)
    else x <- list(...)
    n <- length(x)
    if (n == 0) return(NULL)
    x <- lapply(x, function(y) if (length(y))
        as.matrix(y)
        else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2],
                                                 (y[3] + 1):y[4]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}

.estimated_parameter_names <- function(pmatrix)
{
    estimate <- NULL
    paste0(pmatrix[estimate == 1]$series,".",pmatrix[estimate == 1]$parameter)
}

.check_y_filter <- function(object, y = NULL)
{
    # check index
    if (!is.null(y)) {
        index_new_y <- index(y)
        index_old_y <- object$spec$target$index
        check <- max(index_old_y) < min(index_new_y)
        if (!check) {
            stop("\none of more timestamps in y is before the timestamps in the object data.")
        }
    }
    # check series
    n_series <- object$spec$n_series
    if (n_series != NCOL(y)) {
        stop(paste0("\nexpected ",n_series, " series but got ", NCOL(y)))
    }
    return(y = y)
}


.decorrelate_errors <- function(R, Z)
{
    dims <- dim(R)
    m <- dims[1]
    if (length(dims) == 3) {
        n <- dims[3]
        out <- do.call(rbind, lapply(1:n, function(i){
            e <- eigen(R[,,i], TRUE)
            w <- e$vectors %*% diag(1/sqrt(e$values), m, m) %*% solve(e$vectors)
            matrix(Z[i,] %*% w, nrow = 1)
        }))
    } else {
        e <- eigen(R, TRUE)
        w <- e$vectors %*% diag(1/sqrt(e$values), m, m) %*% solve(e$vectors)
        out <- Z %*% w
    }
    return(out)
}

.check_weights <- function(weights, m, n)
{
    if (is.null(weights)) weights <- matrix(rep(1/m, m), nrow = 1)
    if (!is.matrix(weights)) weights <- matrix(weights, nrow = 1)
    if (NCOL(weights) != m) stop(paste0("\nlength of weights vector must be : ", m))
    if (NROW(weights) > 1) {
        if (NROW(weights) != n) stop(paste0("\nif using a matrix of weights, NROW(weights) must be :", n))
    }
    return(weights)
}

.retain_dimensions_array <- function(x, i)
{
    dims <- dim(x[,,i, drop = FALSE])
    if (is.vector(x[,,i])) {
        return(matrix(x[,,i], nrow = dims[1], ncol = dims[2]))
    } else {
        return(x[,,i])
    }
}

.get_dcc_params <- function(object) {
    group <- NULL
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    shape <- object$parmatrix[group == "shape"]$value
    Qbar <- object$qbar
    Nbar <- object$nbar
    return(list(alpha = alpha, gamma = gamma, beta = beta, shape = shape, Qbar = Qbar, Nbar = Nbar))
}

.get_dcc_order <- function(object) {
    dccorder <- object$spec$dynamics$order
    if (object$spec$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    return(dccorder)
}

.joint_parameter_matrix <- function(spec)
{
    series <- NULL
    out <- lapply(spec$univariate, function(x) x$parmatrix)
    series_names <- names(spec$univariate)
    for (i in 1:length(out)) out[[i]][,series := series_names[i]]
    out <- rbindlist(out)
    out <- rbind(out, spec$parmatrix[,colnames(out), with = FALSE])
    return(out)
}

.get_garch_parameters <- function(object)
{
    parameter <- NULL
    p <- object$parmatrix
    p <- rbind(p[parameter == "skew"]$value,
          p[parameter == "shape"]$value,
          p[parameter == "lambda"]$value)
    rownames(p) <- c("skew","shape","lambda")
    return(p)
}


.rescale_univariate_hessian <- function(x)
{
    D <- diag(1/x$parameter_scale, nrow = length(x$parameter_scale), ncol = length(x$parameter_scale))
    H <- t(D) %*% x$scaled_hessian %*% D
    return(H)
}

.rescale_univariate_scores <- function(x)
{
    estimate <- NULL
    out <- x$scaled_scores
    N <- nrow(out)
    out <- sweep(out, 2, 1/x$parameter_scale, FUN = "*")
    colnames(out) <- x$parmatrix[estimate == 1]$parameter
    return(out)
}

is_square <- function(x)
{
    if (NCOL(x) == NROW(x)) return(TRUE) else return(FALSE)
}

array_matrix_mean <- function(x)
{
    # assumes an s_series x n_series x h x nsim array
    n_series <- dim(x)[2]
    h <- dim(x)[3]
    nsim <- dim(x)[4]
    M <- array(0, dim = c(n_series, n_series, h))
    for (i in 1:h) {
        tmp <- x[,,i,]
        tmp <- apply(tmp, 3, identity, simplify = FALSE)
        M[,,i] <- Reduce(`+`, tmp)/nsim
    }
    return(M)
}

array2distribution <- function(x, pair = c(1,1))
{
    if (length(dim(x)) != 4) stop("\nfunction only valid for 4-d arrays (n_series x n_series x h x nsim")
    indx <- attr(x, "index")
    tmp <- t(x[pair[1], pair[2], , ])
    colnames(tmp) <- indx
    if (is.numeric(indx)) attr(tmp,"date_class") <- "numeric" else attr(tmp,"date_class") <- "Date"
    class(tmp) <- "tsmodel.distribution"
    return(tmp)
}


.check_index_estimate <- function(object, index) {
    index <- as.integer(index)
    max_index <- length(object$spec$target$index)
    if (any(index > max_index)) stop("\nindex > max(timesteps)")
    if (any(index <= 0)) stop("\nindex must be strictly positive.")
    return(index)
}

.check_index_predict <- function(object, index)
{
    index <- as.integer(index)
    max_index <- length(object$univariate[[1]]$sigma)
    if (any(index > max_index)) stop("\nindex > max(horizon)")
    if (any(index <= 0)) stop("\nindex must be strictly positive.")
    return(index)
}

.check_index_simulate <- function(object, index)
{
    index <- as.integer(index)
    max_index <- NCOL(object$univariate[[1]]$sigma)
    if (any(index > max_index)) stop("\nindex > max(horizon)")
    if (any(index <= 0)) stop("\nindex must be strictly positive.")
    return(index)
}

.set_distribution_class <- function(x, index)
{
    colnames(x) <- index
    if (is.character(index)) d_class <- "Date" else d_class <- "numeric"
    attr(x, "date_class") <- d_class
    class(x) <- "tsmodel.distribution"
    return(x)
}

.port_sigma <- function(object, weights)
{
    H <- tscov(object)
    n <- dim(H)[3]
    S <- rep(0, n)
    for (i in 1:n) {
        S[i] <- sqrt(weights[i,] %*% H[,,i] %*% weights[i, ])
    }
    return(S)
}

.zca <- function(x) {
    s <- eigen(x, symmetric = TRUE)
    E <- s$values
    V <- s$vectors
    n <- length(E)
    C <- V %*% diag(1.0/sqrt(E), n, n) %*% t(V)
    return(C)
}

solver_conditions <- function(pars, fn, gr, hess, arglist)
{
    kkttol <- 0.001
    kkt2tol <- 1e-06
    kkt1 <- NA
    kkt2 <- NA
    npar <- length(pars)
    nbm <- 0
    fval <- fn(pars, arglist)
    ngr <- gr(pars, arglist)
    nHes <- hess(pars, arglist)
    pHes <- nHes
    gmax <- max(abs(ngr))
    kkt1 <- (gmax <= kkttol * (1 + abs(fval)))
    phev <- try(eigen(pHes)$values, silent = TRUE)
    if (!inherits(phev, "try-error")) {
        negeig <- (phev[npar] <= (-1) * kkt2tol * (1 + abs(fval)))
        evratio <- phev[npar - nbm]/phev[1]
        kkt2 <- (evratio > kkt2tol) && (!negeig)
        ans <- list(evratio, kkt1, kkt2)
        names(ans) <- c("evratio", "kkt1", "kkt2")
        return(ans)
    }
    else {
        evratio <- NA
        ans <- list(evratio, kkt1, kkt2)
        names(ans) <- c("evratio", "kkt1", "kkt2")
        return(ans)
    }
}

.arma_orders <- function(object)
{
    sapply(object, function(x) {
        arma <- x$spec$model$arma
        if (is.null(arma)) 0L else as.integer(sum(arma))
    })
}

.has_arma_dynamics <- function(object) any(.arma_orders(object) > 0)

.xreg_series <- function(object)
{
    sapply(object, function(x) isTRUE(x$spec$xreg$include_xreg))
}

.has_xreg_dynamics <- function(object) any(.xreg_series(object))

.check_cond_mean_dynamics <- function(object, cond_mean)
{
    if (is.null(cond_mean)) return(invisible(NULL))
    orders <- .arma_orders(object)
    has_arma <- orders > 0
    has_xreg <- .xreg_series(object)
    if (any(has_arma) || any(has_xreg)) {
        snames <- names(object)
        if (is.null(snames)) snames <- paste0("series_", seq_along(orders))
        reasons <- NULL
        if (any(has_arma)) {
            reasons <- c(reasons, paste0("ARMA dynamics (", paste0(snames[has_arma], collapse = ", "), ")"))
        }
        if (any(has_xreg)) {
            reasons <- c(reasons, paste0("mean regressors (", paste0(snames[has_xreg], collapse = ", "), ")"))
        }
        stop(paste0("\ncond_mean cannot be used when any first stage model has ",
                    paste0(reasons, collapse = " or "),
                    ". Either use arma = c(0,0) with no mean regressors in all first stage models and supply cond_mean, or let the first stage models generate the conditional mean."), call. = FALSE)
    }
    return(invisible(NULL))
}

# validates and normalises the mean equation regressor argument used by
# tsfilter (newxreg), predict (newxreg) and simulate (xreg). The user supplies
# a list with one element per series; each element is a matrix/xts of
# regressors or NULL. Returns a list of length n_series whose elements are
# either NULL or a numeric matrix with n_required rows.
.multi_xreg_spec <- function(univariate, newxreg, n_required, argument_name = "newxreg")
{
    n_series <- length(univariate)
    snames <- names(univariate)
    if (is.null(snames)) snames <- paste0("series_", seq_len(n_series))
    has_xreg <- .xreg_series(univariate)
    ncols <- sapply(seq_len(n_series), function(i) {
        if (has_xreg[i]) NCOL(univariate[[i]]$spec$xreg$xreg) else 0L
    })
    zero_fill <- function(i) matrix(0, nrow = n_required, ncol = ncols[i])
    zero_filled <- character(0)
    ignored <- character(0)
    if (is.null(newxreg)) {
        if (!any(has_xreg)) return(vector("list", n_series))
        zero_filled <- snames[has_xreg]
        out <- lapply(seq_len(n_series), function(i) if (has_xreg[i]) zero_fill(i) else NULL)
        warning(paste0("\n", argument_name, " is NULL but the following series were specified with mean regressors: ",
                       paste0(zero_filled, collapse = ", "), ". Setting to zero."), call. = FALSE)
        return(out)
    }
    if (!is.list(newxreg)) {
        stop(paste0("\n", argument_name, " must be a list with one element per series (a matrix/xts of regressors or NULL)."), call. = FALSE)
    }
    nms <- names(newxreg)
    if (is.null(nms) || all(!nzchar(nms))) {
        if (length(newxreg) != n_series) {
            stop(paste0("\nan unnamed ", argument_name, " list must have exactly one element per series (", n_series, " expected, got ", length(newxreg), ")."), call. = FALSE)
        }
        xreg_list <- newxreg
    } else if (all(nzchar(nms))) {
        if (any(duplicated(nms))) {
            stop(paste0("\n", argument_name, " has duplicated names: ", paste0(unique(nms[duplicated(nms)]), collapse = ", ")), call. = FALSE)
        }
        unknown <- setdiff(nms, snames)
        if (length(unknown) > 0) {
            stop(paste0("\n", argument_name, " contains names not matching any first stage series: ",
                        paste0(unknown, collapse = ", ")), call. = FALSE)
        }
        xreg_list <- vector("list", n_series)
        xreg_list[match(nms, snames)] <- newxreg
    } else {
        stop(paste0("\n", argument_name, " must be either a fully named list (one element per series name) or a fully unnamed list (one element per series)."), call. = FALSE)
    }
    out <- vector("list", n_series)
    for (i in seq_len(n_series)) {
        x <- xreg_list[[i]]
        if (is.null(x)) {
            if (has_xreg[i]) {
                out[[i]] <- zero_fill(i)
                zero_filled <- c(zero_filled, snames[i])
            }
        } else {
            if (!has_xreg[i]) {
                ignored <- c(ignored, snames[i])
            } else {
                if (is.xts(x)) x <- coredata(x)
                x <- as.matrix(x)
                if (NROW(x) != n_required) {
                    stop(paste0("\n", argument_name, " for series ", snames[i], " must have ", n_required, " rows (got ", NROW(x), ")."), call. = FALSE)
                }
                if (NCOL(x) != ncols[i]) {
                    stop(paste0("\n", argument_name, " for series ", snames[i], " must have ", ncols[i], " columns (got ", NCOL(x), ")."), call. = FALSE)
                }
                if (any(!is.finite(x))) {
                    stop(paste0("\nNA/NaN/Inf values found in ", argument_name, " for series ", snames[i], "."), call. = FALSE)
                }
                out[[i]] <- x
            }
        }
    }
    if (length(zero_filled) > 0) {
        warning(paste0("\nno ", argument_name, " supplied for the following series which were specified with mean regressors: ",
                       paste0(zero_filled, collapse = ", "), ". Setting to zero."), call. = FALSE)
    }
    if (length(ignored) > 0) {
        warning(paste0("\n", argument_name, " was supplied for the following series which were not specified with mean regressors: ",
                       paste0(ignored, collapse = ", "), ". Ignoring."), call. = FALSE)
    }
    return(out)
}

.cond_mean_spec <- function(mu = NULL, n_series, n_points, series_names)
{
    if (is.null(mu)) {
        mu <- matrix(0, ncol = n_series, nrow = n_points)
        colnames(mu) <- series_names
    } else {
        if (!is.matrix(mu)) stop(paste0("\ncond_mean must be a matrix of dimenions ", n_points, "x ", n_series))
        if (NROW(mu) != n_points) stop(paste0("\ncond_mean must be a matrix of dimenions ", n_points, "x ", n_series))
        if (NCOL(mu) != n_series) stop(paste0("\ncond_mean must be a matrix of dimenions ", n_points, "x ", n_series))
        mu <- coredata(mu)
        colnames(mu) <- series_names
    }
    return(mu)
}


.cond_mean_inject <- function(x, value, recenter = FALSE) {
    n_series <- dim(x)[2]
    h <- dim(x)[1]
    for (i in 1:n_series) {
        if (recenter) {
            if (h == 1) {
                x[,i,] <- x[,i,] - mean(x[,i,])
            } else {
                x[,i,] <- sweep(x[,i,], 1, rowMeans(x[,i,]), FUN = "-")
            }
        }
        if (h == 1) {
            x[,i,] <- x[,i,] + value[,i]
        } else {
            x[,i,] <- sweep(x[,i,], 1, value[,i], FUN = "+")
        }
    }
    return(x)
}

.covariance_prediction_plots <- function(model, p, series = c(1,2))
{
    hist <- tscov(model)[series[1], series[2],]
    hist <- xts(hist, order.by = model$spec$target$index)
    pred <- t(tscov(p)[series[1], series[2],,])
    if (grepl("predict",class(p)[1])) {
        colnames(pred) <- p$future_dates
        class(pred) <- "tsmodel.distribution"
    } else {
        nh <- as.integer(median(diff(index(hist))))
        f_dates <- as.Date(tail(model$spec$target$index, 1)) + seq(nh, nh * NCOL(pred), by = nh)
        colnames(pred) <- as.character(f_dates)
        class(pred) <- "tsmodel.distribution"
    }
    L <- list(distribution = pred, original_series = hist)
    class(L) <- "tsmodel.predict"
    return(L)
}


.es_empirical_calculation <- function(x, alpha = 0.05) {
    x <- sort(x)
    var_threshold <- quantile(x, alpha)
    expected_shortfall <- mean(x[x < var_threshold])
    return(expected_shortfall)
}

.make_xts <- function(x, dates)
{
   return(xts(x, order.by = dates))
}
