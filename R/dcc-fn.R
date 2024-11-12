# !diagnostics suppress=data.table,copy,as.data.table,setnames

# dcc c++ values ---------------------------------------------------

.dcc_dynamic_values <- function(pars, spec, type = "nll", return_all = FALSE)
{
    estimate <- group <- NULL
    type <- match.arg(type[1], c("nll","dcc_nll","ll_vec","dccll_vec","R","Qbar","Nbar","Q"))
    pmatrix <- copy(spec$parmatrix)
    pmatrix[estimate == 1]$value <- pars
    dist <- spec$distribution
    dccorder <- spec$dynamics$order
    if (spec$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    z <- residuals(spec$univariate, standardize = TRUE)
    s <- coredata(sigma(spec$univariate))
    if (return_all) {
        out <- switch(dist,
                      "mvn" = .dcc_dynamic_normal(alpha = pmatrix[group == "alpha"]$value,
                                                       gamma = pmatrix[group == "gamma"]$value,
                                                       beta = pmatrix[group == "beta"]$value,
                                                       z = z, s = s,
                                                       dccorder = dccorder),
                      "mvt" = .dcc_dynamic_student(alpha = pmatrix[group == "alpha"]$value,
                                                       gamma = pmatrix[group == "gamma"]$value,
                                                       beta = pmatrix[group == "beta"]$value,
                                                       shape = pmatrix[group == "shape"]$value,
                                                       z = z, s = s,
                                                       dccorder = dccorder)
        )
    } else {
        out <- switch(dist,
                      "mvn" = .dcc_dynamic_normal(alpha = pmatrix[group == "alpha"]$value,
                                                       gamma = pmatrix[group == "gamma"]$value,
                                                       beta = pmatrix[group == "beta"]$value,
                                                       z = z, s = s,
                                                       dccorder = dccorder)[[type]],
                      "mvt" = .dcc_dynamic_student(alpha = pmatrix[group == "alpha"]$value,
                                                       gamma = pmatrix[group == "gamma"]$value,
                                                       beta = pmatrix[group == "beta"]$value,
                                                       shape = pmatrix[group == "shape"]$value,
                                                       z = z, s = s,
                                                       dccorder = dccorder)[[type]]
        )
    }
    return(out)
}

.dcc_constant_values <- function(pars, spec, type = "nll", return_all = FALSE)
{
    estimate <- group <- NULL
    pmatrix <- copy(spec$parmatrix)
    if (!is.null(pars)) {
        pmatrix[estimate == 1]$value <- pars
    }
    dist <- spec$distribution
    z <- residuals(spec$univariate, standardize = TRUE)
    s <- coredata(sigma(spec$univariate))
    if (return_all) {
        out <- switch(dist,
                      "mvn" = .dcc_constant_normal(Z = z, S = s),
                      "mvt" = .dcc_constant_student(Z = z, S = s, shape = pmatrix[group == "shape"]$value)
        )
    } else {
        out <- switch(dist,
                      "mvn" = .dcc_constant_normal(Z = z, S = s)[[type]],
                      "mvt" = .dcc_constant_student(Z = z, S = s, shape = pmatrix[group == "shape"]$value)[[type]]
        )
    }
    return(out)
}

# dcc c++ fun ---------------------------------------------------
.dcc_dynamic_fun <- function(spec, type = "nll")
{
    estimate <- group <- NULL
    if (spec$distribution == "mvn") {
        fun <- function(x, arglist)
        {
            arglist$parmatrix[estimate == 1]$value <- x
            .dcc_dynamic_normal(alpha = arglist$parmatrix[group == "alpha"]$value,
                                gamma = arglist$parmatrix[group == "gamma"]$value,
                                beta = arglist$parmatrix[group == "beta"]$value,
                                z = arglist[["z"]], s = arglist[["s"]], dccorder = arglist$dccorder)[[type]]
        }
    } else {
        fun <- function(x, arglist)
        {
            arglist$parmatrix[estimate == 1]$value <- x
            .dcc_dynamic_student(alpha = arglist$parmatrix[group == "alpha"]$value,
                                 gamma = arglist$parmatrix[group == "gamma"]$value,
                                 beta = arglist$parmatrix[group == "beta"]$value,
                                 shape = arglist$parmatrix[group == "shape"]$value,
                                 z = arglist[["z"]], s = arglist[["s"]], dccorder = arglist$dccorder)[[type]]
        }
    }
    if (spec$dynamics$model == "adcc") {
        cons <- function(x, arglist) {
            arglist$parmatrix[estimate == 1]$value <- x
            .adcc_constraint(alpha = arglist$parmatrix[group == "alpha"]$value,
                             gamma = arglist$parmatrix[group == "gamma"]$value,
                             beta = arglist$parmatrix[group == "beta"]$value,
                             shape = arglist$parmatrix[group == "shape"]$value,
                             z = arglist[["z"]],
                             dccorder = arglist$dccorder) - 0.999
        }
        jac <- function(x, arglist) {
            jacobian(func = arglist$inequality_cons, x = x, arglist = arglist)
        }
    } else {
        cons <- function(x, arglist) {
            arglist$parmatrix[estimate == 1]$value <- x
            sum(arglist$parmatrix[group == "alpha"]$value) + sum(arglist$parmatrix[group == "beta"]$value) - 0.999
        }

        j <- matrix(0, ncol = length(spec$parmatrix[estimate == 1]$value), nrow = 1)
        cnames <- spec$parmatrix[estimate == 1]$group
        if (any(cnames %in% "alpha")) {
            j[which(cnames %in% "alpha")] <- 1
        }
        if (any(cnames %in% "beta")) {
            j[which(cnames %in% "beta")] <- 1
        }
        jac <- function(x, arglist) {
            return(j)
        }
    }
    grad_fun <- function(x, arglist) {
        if (arglist$inequality_cons(x, arglist) >= 0) {
            return(matrix(1e6, ncol = length(x), nrow = 1))
        } else {
            return(grad(func = arglist$fun, x = x, arglist = arglist))
        }
    }
    hess_fun <- function(x, arglist) {
        nderiv_hessian(func = arglist$fun, x = x, lower = arglist$lower, upper = arglist$upper, arglist = arglist)
    }
    return(list(fun = fun, grad = grad_fun, hess = hess_fun, inequality_cons = cons, inequality_jac = jac))
}


.dcc_constant_fun <- function(spec, type = "nll")
{
    estimate <- group <- NULL
    if (spec$distribution == "mvt") {
        fun <- function(x, arglist)
        {
            arglist$parmatrix[estimate == 1]$value <- x
            .dcc_constant_student(shape = arglist$parmatrix[group == "shape"]$value, S = arglist[["s"]], Z = arglist$z)[[type]]
        }
    } else {
        stop("\ndistrbution not recognized.")
    }
    grad_fun <- function(x, arglist) {
        return(grad(func = arglist$fun, x = x, arglist = arglist))
    }
    hess_fun <- function(x, arglist) {
        nderiv_hessian(func = arglist$fun, x = x, lower = arglist$lower, upper = arglist$upper, arglist = arglist)
    }
    return(list(fun = fun, grad = grad_fun, hess = hess_fun, inequality_cons = NULL, inequality_jac = NULL))
}

.dcc_dynamic_joint_nll <- function(pars, spec)
{
    estimate <- NULL
    pmatrix <- spec$joint_parmatrix
    pmatrix[estimate == 1]$value <- pars
    pmatrix_spec <- split(pmatrix, by = "series")
    pmatrix_dcc <- pmatrix_spec[["DCC"]]
    pmatrix_spec[["DCC"]] <- NULL
    n <- length(pmatrix_spec)
    new_fit <- lapply(1:n, function(i){
        newf <- spec$univariate[[i]]
        newf$parmatrix <- pmatrix_spec[[i]]
        maxpq <- max(newf$spec$model$order)
        newf$sigma <- spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$sigma
        if (maxpq > 0) newf$sigma <- newf$sigma[-c(1:maxpq)]
        newf$loglik <- spec$univariate[[i]]$tmb$fn(pmatrix_spec[[i]][estimate == 1]$value)
        return(newf)
    })
    new_fit <- to_multi_estimate(new_fit)
    spec$univariate <- new_fit
    tmp <- .dcc_dynamic_values(pmatrix_dcc[estimate == 1]$value, spec, type = "nll")
    model_loglik <- -1.0 * tmp
    return(model_loglik)
}

.dcc_dynamic_joint_nllvec <- function(pars, spec)
{
    estimate <- NULL
    pmatrix <- spec$joint_parmatrix
    pmatrix[estimate == 1]$value <- pars
    pmatrix_spec <- split(pmatrix, by = "series")
    pmatrix_dcc <- pmatrix_spec[["DCC"]]
    pmatrix_spec[["DCC"]] <- NULL
    n <- length(pmatrix_spec)
    new_fit <- lapply(1:n, function(i) {
        newf <- spec$univariate[[i]]
        newf$parmatrix <- pmatrix_spec[[i]]
        maxpq <- max(newf$spec$model$order)
        newf$sigma <- spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$sigma
        if (maxpq > 0) newf$sigma <- newf$sigma[-c(1:maxpq)]
        newf$llvec <- -1 * log(spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$ll_vector)
        return(newf)
    })
    new_fit <- to_multi_estimate(new_fit)
    spec$univariate <- new_fit
    model_loglik <- -1.0 * .dcc_dynamic_values(pmatrix_dcc[estimate == 1]$value, spec, type = "dccll_vec")
    maxpq <- max(spec$dynamics$order)
    if (maxpq > 0) {
        model_loglik <- model_loglik[-c(1:maxpq)]
    }
    return(model_loglik)
}


.dcc_dynamic_hessian <- function(pars, spec)
{
    estimate <- NULL
    garch_hessian <- lapply(spec$univariate, function(x) .rescale_univariate_hessian(x))
    m <- length(spec$parmatrix[estimate == 1]$value)
    dcc_hessian <- matrix(0, m, m)
    H <- bdiag(garch_hessian)
    H <- bdiag(H, dcc_hessian)
    init_llh <- .dcc_dynamic_joint_nll(pars, spec)
    .epsilon <- .Machine$double.eps
    n <- length(pars)
    step_size <- .epsilon^(1/3) * pmax(abs(pars), 1)
    tmp <- pars + step_size
    step_size <- tmp - pars
    deltas <- as.matrix(diag(step_size, length(step_size), length(step_size)))
    g <- NULL
    g <- future_lapply(1:n, function(i) {
        setDTthreads(1)
        .dcc_dynamic_joint_nll(pars + deltas[, i], spec)
    }, future.packages = c("tsmarch","data.table"), future.seed = TRUE)
    g <- eval(g)
    g <- unlist(g)
    step_matrix <- step_size %*% t(step_size)
    tmp <- NULL
    tmp <- future_lapply(1:n, function(i) {
        s_tmp = step_matrix
        for (j in (n - m + 1):n) {
            if (i <= j) {
                s_tmp[i, j] <- (.dcc_dynamic_joint_nll(pars + deltas[, i] + deltas[, j], spec) - g[i] - g[j] + init_llh) / s_tmp[i, j]
                s_tmp[j, i] <- s_tmp[i, j]
            }
        }
        return(s_tmp)
    }, future.packages = c("tsmarch"), future.seed = TRUE)
    tmp <- eval(tmp)
    for (i in 1:n) {
        for (j in (n - m + 1):n) {
            if (i <= j) {
                step_matrix[i, j] = tmp[[i]][i, j]
                step_matrix[j, i] = step_matrix[i, j]
            }
        }
    }
    lower_block <- step_matrix[(n - m + 1):n, ]
    H[(n - m + 1):n, ] = lower_block
    return(H)
}

.dcc_dynamic_scores <- function(pars, spec)
{
    estimate <- NULL
    garch_scores <- do.call(cbind, lapply(spec$univariate, function(x) .rescale_univariate_scores(x)))
    N <- NROW(garch_scores)
    n <- length(pars)
    m <- length(spec$parmatrix[estimate == 1]$value)
    .epsilon <- .Machine$double.eps
    dcc_pars <- tail(pars, m)
    step_size <-  pmax(abs(dcc_pars/2), 0.01) * .epsilon^(1/3)
    step_plus <- dcc_pars + step_size
    step_minus <- dcc_pars - step_size
    likvec_plus <- likvec_minus <- matrix(0, ncol = m, nrow = N)
    pars_plus <- pars_minus <- pars
    for (i in 1:m) {
        step_pars_plus <- dcc_pars
        step_pars_minus <- dcc_pars
        step_pars_plus[i] <- step_plus[i]
        step_pars_minus[i] <- step_minus[i]
        pars_plus[(n - m + 1):n] <- step_pars_plus
        pars_minus[(n - m + 1):n] <- step_pars_minus
        likvec_plus[,i] <- .dcc_dynamic_joint_nllvec(pars_plus, spec)
        likvec_minus[,i] <- .dcc_dynamic_joint_nllvec(pars_minus, spec)
    }
    central_diff <- likvec_plus - likvec_minus
    dcc_scores <- matrix(NA, ncol = m, nrow = N)
    central_diff_sd = 2 * kronecker(matrix(1, N, 1), t(step_size))
    for (i in 1:m) {
        dcc_scores[, i] = central_diff[, i]/central_diff_sd[, i]
    }
    joint_scores <- cbind(garch_scores, dcc_scores)
    return(joint_scores)
}


.dcc_constant_joint_nll <- function(pars, spec)
{
    estimate <- NULL
    pmatrix <- spec$joint_parmatrix
    pmatrix[estimate == 1]$value <- pars
    pmatrix_spec <- split(pmatrix, by = "series")
    pmatrix_dcc <- pmatrix_spec[["CONSTANT"]]
    pmatrix_spec[["CONSTANT"]] <- NULL
    n <- length(pmatrix_spec)
    new_fit <- lapply(1:n, function(i){
        newf <- spec$univariate[[i]]
        newf$parmatrix <- pmatrix_spec[[i]]
        maxpq <- max(newf$spec$model$order)
        newf$sigma <- spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$sigma
        if (maxpq > 0) newf$sigma <- newf$sigma[-c(1:maxpq)]
        newf$loglik <- spec$univariate[[i]]$tmb$fn(pmatrix_spec[[i]][estimate == 1]$value)
        return(newf)
    })
    new_fit <- to_multi_estimate(new_fit)
    spec$univariate <- new_fit
    model_loglik <- .dcc_constant_values(pmatrix_dcc[estimate == 1]$value, spec, type = "nll")
    model_loglik <- -1.0 * model_loglik
    return(model_loglik)
}


.dcc_constant_joint_nllvec <- function(pars, spec)
{
    estimate <- NULL
    pmatrix <- spec$joint_parmatrix
    pmatrix[estimate == 1]$value <- pars
    pmatrix_spec <- split(pmatrix, by = "series")
    pmatrix_dcc <- pmatrix_spec[["CONSTANT"]]
    pmatrix_spec[["CONSTANT"]] <- NULL
    n <- length(pmatrix_spec)
    new_fit <- lapply(1:n, function(i){
        newf <- spec$univariate[[i]]
        newf$parmatrix <- pmatrix_spec[[i]]
        maxpq <- max(newf$spec$model$order)
        newf$sigma <- spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$sigma
        if (maxpq > 0) newf$sigma <- newf$sigma[-c(1:maxpq)]
        newf$llvec <- -1 * log(spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$ll_vector)
        return(newf)
    })
    new_fit <- to_multi_estimate(new_fit)
    spec$univariate <- new_fit
    model_loglik <- -1.0 * .dcc_constant_values(pmatrix_dcc[estimate == 1]$value, spec, type = "ll_vec")
    return(model_loglik)
}

.dcc_constant_hessian <- function(pars, spec)
{
    estimate <- NULL
    garch_hessian <- lapply(spec$univariate, function(x) .rescale_univariate_hessian(x))
    m <- length(spec$parmatrix[estimate == 1]$value)
    dcc_hessian <- matrix(0, m, m)
    H <- bdiag(garch_hessian)
    H <- bdiag(H, dcc_hessian)
    init_llh <- .dcc_constant_joint_nll(pars, spec)
    .epsilon <- .Machine$double.eps
    n <- length(pars)
    step_size <- .epsilon^(1/3) * pmax(abs(pars), 1)
    tmp <- pars + step_size
    step_size <- tmp - pars
    deltas <- as.matrix(diag(step_size, length(step_size), length(step_size)))
    g <- NULL
    g <- future_lapply(1:n, function(i) {
        setDTthreads(1)
        .dcc_constant_joint_nll(pars + deltas[, i], spec)
    }, future.packages = c("tsmarch","data.table"), future.seed = NULL)
    g <- eval(g)
    g <- unlist(g)
    step_matrix <- step_size %*% t(step_size)
    tmp <- NULL
    tmp <- future_lapply(1:n, function(i) {
        s_tmp = step_matrix
        for (j in (n - m + 1):n) {
            if (i <= j) {
                s_tmp[i, j] <- (.dcc_constant_joint_nll(pars + deltas[, i] + deltas[, j], spec) - g[i] - g[j] + init_llh) / s_tmp[i, j]
                s_tmp[j, i] <- s_tmp[i, j]
            }
        }
        return(s_tmp)
    }, future.packages = c("tsmarch"), future.seed = TRUE)
    tmp <- eval(tmp)
    for (i in 1:n) {
        for (j in (n - m + 1):n) {
            if (i <= j) {
                step_matrix[i, j] = tmp[[i]][i, j]
                step_matrix[j, i] = step_matrix[i, j]
            }
        }
    }
    lower_block <- step_matrix[(n - m + 1):n, ]
    H[(n - m + 1):n, ] = lower_block
    return(H)
}

.dcc_constant_scores <- function(pars, spec)
{
    estimate <- NULL
    garch_scores <- do.call(cbind, lapply(spec$univariate, function(x) .rescale_univariate_scores(x)))
    N <- NROW(garch_scores)
    n <- length(pars)
    m <- length(spec$parmatrix[estimate == 1]$value)
    .epsilon <- .Machine$double.eps
    dcc_pars <- tail(pars, m)
    step_size <-  pmax(abs(dcc_pars/2), 0.01) * .epsilon^(1/3)
    step_plus <- dcc_pars + step_size
    step_minus <- dcc_pars - step_size
    likvec_plus <- likvec_minus <- matrix(0, ncol = m, nrow = N)
    pars_plus <- pars_minus <- pars
    for (i in 1:m) {
        step_pars_plus <- dcc_pars
        step_pars_minus <- dcc_pars
        step_pars_plus[i] <- step_plus[i]
        step_pars_minus[i] <- step_minus[i]
        pars_plus[(n - m + 1):n] <- step_pars_plus
        pars_minus[(n - m + 1):n] <- step_pars_minus
        likvec_plus[,i] <- .dcc_constant_joint_nllvec(pars_plus, spec)
        likvec_minus[,i] <- .dcc_constant_joint_nllvec(pars_minus, spec)
    }
    central_diff <- likvec_plus - likvec_minus
    dcc_scores <- matrix(NA, ncol = m, nrow = N)
    central_diff_sd = 2 * kronecker(matrix(1, N, 1), t(step_size))
    for (i in 1:m) {
        dcc_scores[, i] = central_diff[, i]/central_diff_sd[, i]
    }
    joint_scores <- cbind(garch_scores, dcc_scores)
    return(joint_scores)
}

.dcc_sim_dynamic_initialize_values <- function(object, Q_init, Z_init, init_method)
{
    maxpq <- max(object$spec$dynamics$order)
    n_series <- object$spec$n_series
    Qinit <- array(0, dim = c(n_series, n_series, min(1,maxpq)))
    Zinit <- matrix(0, ncol = n_series, nrow = min(1,maxpq))
    if (maxpq > 0) {
        if (!is.null(Q_init)) {
            if (!all.equal(dim(Q_init), dim(Qinit))) stop("\nQ_init must be an array of size (n_series, n_series, maxpq)")
            Qinit <- Q_init
        } else {
            if (init_method == "start") {
                for (i in 1:maxpq) {
                    Qinit[,,i] <- object$Qbar
                }
            } else if (init_method == "end") {
                N <- NROW(object$Q)
                dims <- c(n_series, n_series, N)
                Q <- .compose_tri_matrix(object$Q, dims = dims, diag = TRUE)
                for (i in 1:maxpq) {
                    Qinit[,,i] <- Q[,,(N - maxpq + i)]
                }
            }
        }
        if (!is.null(Z_init)) {
            if (!all.equal(dim(Z_init), dim(Zinit))) stop("\nZ_init must be a matrix of size (maxpq, n_series)")
            Zinit <- Z_init
        } else {
            Z <- residuals(object, standardize = TRUE)
            if (init_method == "start") {
                for (i in 1:maxpq) {
                    Zinit[i, ] <- matrix(0, ncol = n_series, nrow = maxpq)
                }
            } else if (init_method == "end") {
                N <- NROW(Z)
                for (i in 1:maxpq) {
                    Zinit[i,] <- Z[(N - maxpq + i), ]
                }
            }
        }
    }
    return(list(Qinit = Qinit, Zinit = Zinit))
}
