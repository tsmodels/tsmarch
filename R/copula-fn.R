# !diagnostics suppress=data.table,copy,as.data.table,setnames

.copula_parameters <- function(dynamics, copula, order = c(1,1))
{
    include <- series <- out  <- NULL
    if (dynamics == "constant") {
        model_name <- "CONSTANT"
        if (copula == "mvt") {
            out <- data.table(parameter = "shape", value = 4.0, lower = 2.01, upper = 50, estimate = 1, scale = 1, group = "shape", equation = "[D]", symbol = paste0("\\nu"))
        } else {
            out <- data.table(parameter = "shape", value = 4.0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "shape", equation = "[D]", symbol = paste0("\\nu"))
        }
        out[,include := 1]
        out[,series := model_name]
    } else {
        model_name <- "DCC"
        if (order[1] > 0) {
            out <- data.table(parameter = paste0("alpha_",1:order[1]), value = 0.05, lower = 0, upper = 1, estimate = ifelse(dynamics == "constant", 0, 1), scale = 1, group = "alpha", equation = "[DCC]", symbol = paste0("\\alpha_",1:order[1]))
            if (dynamics == "adcc") {
                tmp <- data.table(parameter = paste0("gamma_",1:order[1]), value = 0.02, lower = 0, upper = 1, estimate = 1, scale = 1, group = "gamma", equation = "[DCC]", symbol = paste0("\\gamma_",1:order[1]))
            } else {
                tmp <- data.table(parameter = paste0("gamma_",1:order[1]), value = 0.0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "gamma", equation = "[DCC]", symbol = paste0("\\gamma_",1:order[1]))
            }
            out <- rbind(out, tmp)
        } else {
            out <- data.table(parameter = paste0("alpha_",1), value = 0.0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "alpha", equation = "[DCC]", symbol = paste0("\\alpha_",1))
        }
        if (order[2] > 0) {
            tmp <- data.table(parameter = paste0("beta_",1:order[1]), value = 0.8, lower = 0, upper = 1, estimate = ifelse(dynamics == "constant", 0, 1), scale = 1, group = "beta",  equation = "[DCC]", symbol = paste0("\\beta_",1:order[2]))
        } else{
            tmp <- data.table(parameter = paste0("beta_",1), value = 0.0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "beta", equation = "[DCC]", symbol = paste0("\\alpha_",1))
        }
        out <- rbind(out, tmp)
        if (copula == "mvt") {
            tmp <- data.table(parameter = "shape", value = 4.0, lower = 2.01, upper = 50, estimate = 1, scale = 1, group = "shape", equation = "[D]", symbol = paste0("\\nu"))
        } else {
            tmp <- data.table(parameter = "shape", value = 4.0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "shape", equation = "[D]", symbol = paste0("\\nu"))
        }
        out <- rbind(out, tmp)
    }
    out[,include := 1]
    out[,series := model_name]
    return(out)
}

.copula_ptransform <- function(z, shape, distribution)
{
    if (distribution == "mvn") {
        u <- do.call(cbind, lapply(1:ncol(z), function(i) pnorm(z[,i])))
    } else if (distribution == "mvt") {
        u <- do.call(cbind, lapply(1:ncol(z), function(i) pdist("std", z[,i], shape = shape)))
    }
    return(u)
}

.copula_qtransform <- function(object, u, k)
{
    parameter <- NULL
    transformation <- object$spec$transformation
    if (transformation == "parametric") {
        s <- object$spec$univariate[[k]]
        parameter <- NULL
        dist <- s$spec$distribution
        skew <- s$parmatrix[parameter == "skew"]$value
        shape <- s$parmatrix[parameter == "shape"]$value
        lambda <- s$parmatrix[parameter == "lambda"]$value
        out <- do.call(cbind, lapply(1:NCOL(u), function(i) qdist(distribution = dist, u[,i,drop = FALSE], mu = 0, sigma = 1, skew = skew, shape = shape, lambda = lambda)))
    } else if (transformation == "empirical") {
        ecd_model <- object$spec$transform$transform_model[[k]]
        out <- do.call(cbind, lapply(1:NCOL(u), function(i) quantile(ecd_model, u[,i])))
    } else {
        spd_model <- object$spec$transform$transform_model[[k]]
        out <- do.call(cbind, lapply(1:NCOL(u), function(i) qspd(p = u[,i,drop = FALSE], object = spd_model)))
    }
    return(out)
}


.copula_dynamic_joint_nll <- function(pars, spec)
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
    tmp <- copula_transformation_estimate(new_fit, spec$transformation)
    spec$transform$u <- tmp$u
    spec$transform$transform_model <- tmp$transform_model
    garch_loglik <- sapply(new_fit, function(x) x$loglik)
    copula_loglik <- .copula_dynamic_values(pmatrix_dcc[estimate == 1]$value, spec, type = "nll")
    model_loglik <- -1.0 * (sum(garch_loglik) + copula_loglik)
    return(model_loglik)
}

.copula_dynamic_joint_nllvec <- function(pars, spec)
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
        newf$llvec <- -1 * log(spec$univariate[[i]]$tmb$report(pmatrix_spec[[i]][estimate == 1]$value)$ll_vector)
        return(newf)
    })
    new_fit <- to_multi_estimate(new_fit)
    tmp <- copula_transformation_estimate(new_fit, spec$transformation)
    spec$transform$u <- tmp$u
    spec$transform$transform_model <- tmp$transform_model
    spec$transform$u[spec$transform$u < 3.330669e-16] <- 2.220446e-16
    spec$transform$u[spec$transform$u > 0.99999] <- 0.99999
    copula_loglik <- .copula_dynamic_values(pmatrix_dcc[estimate == 1]$value, spec, type = "ll_vec")
    maxpq <- max(spec$dynamics$order)
    if (maxpq > 0) {
        copula_loglik <- copula_loglik[-c(1:maxpq)]
    }
    garch_loglik <- do.call(cbind, lapply(new_fit, function(x) x$llvec))
    model_loglik <- -1.0 * rowSums(cbind(garch_loglik, copula_loglik))
    return(model_loglik)
}


.copula_dynamic_hessian <- function(pars, spec)
{
    estimate <- NULL
    garch_hessian <- lapply(spec$univariate, function(x) .rescale_univariate_hessian(x))
    m <- length(spec$parmatrix[estimate == 1]$value)
    dcc_hessian <- matrix(0, m, m)
    H <- bdiag(garch_hessian)
    H <- bdiag(H, dcc_hessian)
    init_llh <- .copula_dynamic_joint_nll(pars, spec)
    .epsilon <- .Machine$double.eps
    n <- length(pars)
    step_size <- .epsilon^(1/3) * pmax(abs(pars), 1)
    tmp <- pars + step_size
    step_size <- tmp - pars
    deltas <- as.matrix(diag(step_size, length(step_size), length(step_size)))
    g <- NULL
    g <- future_lapply(1:n, function(i) {
        setDTthreads(1)
        .copula_dynamic_joint_nll(pars + deltas[, i], spec)
    }, future.packages = c("tsmarch"), future.seed = TRUE)
    g <- eval(g)
    g <- unlist(g)
    step_matrix <- step_size %*% t(step_size)
    tmp <- NULL
    tmp <- future_lapply(1:n, function(i) {
        s_tmp = step_matrix
        for (j in (n - m + 1):n) {
            if (i <= j) {
                s_tmp[i, j] <- (.copula_dynamic_joint_nll(pars + deltas[, i] + deltas[, j], spec) - g[i] - g[j] + init_llh) / s_tmp[i, j]
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

.copula_dynamic_scores <- function(pars, spec)
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
        likvec_plus[,i] <- .copula_dynamic_joint_nllvec(pars_plus, spec)
        likvec_minus[,i] <- .copula_dynamic_joint_nllvec(pars_minus, spec)
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


.copula_constant_joint_nll <- function(pars, spec)
{
    estimate <- NULL
    pmatrix <- spec$joint_parmatrix
    pmatrix[estimate == 1]$value <- pars
    pmatrix_spec <- split(pmatrix, by = "series")
    pmatrix_copula <- pmatrix_spec[["CONSTANT"]]
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
    tmp <- copula_transformation_estimate(new_fit, spec$transformation)
    spec$transform$u <- tmp$u
    spec$transform$transform_model <- tmp$transform_model
    garch_loglik <- sapply(new_fit, function(x) x$loglik)
    copula_loglik <- .copula_constant_values(pmatrix_copula[estimate == 1]$value, spec, type = "nll")
    model_loglik <- -1.0 * (sum(garch_loglik) + copula_loglik)
    return(model_loglik)
}


.copula_constant_joint_nllvec <- function(pars, spec)
{
    estimate <- NULL
    pmatrix <- spec$joint_parmatrix
    pmatrix[estimate == 1]$value <- pars
    pmatrix_spec <- split(pmatrix, by = "series")
    pmatrix_copula <- pmatrix_spec[["CONSTANT"]]
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
    tmp <- copula_transformation_estimate(new_fit, spec$transformation)
    spec$transform$u <- tmp$u
    spec$transform$transform_model <- tmp$transform_model
    copula_loglik <- .copula_constant_values(pmatrix_copula[estimate == 1]$value, spec, type = "ll_vec")
    garch_loglik <- do.call(cbind, lapply(new_fit, function(x) x$llvec))
    model_loglik <- -1.0 * rowSums(cbind(garch_loglik, copula_loglik))
    return(model_loglik)
}

.copula_constant_hessian <- function(pars, spec)
{
    estimate <- NULL
    garch_hessian <- lapply(spec$univariate, function(x) .rescale_univariate_hessian(x))
    m <- length(spec$parmatrix[estimate == 1]$value)
    copula_hessian <- matrix(0, m, m)
    H <- bdiag(garch_hessian)
    H <- bdiag(H, copula_hessian)
    init_llh <- .copula_constant_joint_nll(pars, spec)
    .epsilon <- .Machine$double.eps
    n <- length(pars)
    step_size <- .epsilon^(1/3) * pmax(abs(pars), 1)
    tmp <- pars + step_size
    step_size <- tmp - pars
    deltas <- as.matrix(diag(step_size, length(step_size), length(step_size)))
    g <- NULL
    g <- future_lapply(1:n, function(i) {
        setDTthreads(1)
        .copula_constant_joint_nll(pars + deltas[, i], spec)
    }, future.packages = c("tsmarch"), future.seed = NULL)
    g <- eval(g)
    g <- unlist(g)
    step_matrix <- step_size %*% t(step_size)
    tmp <- NULL
    tmp <- future_lapply(1:n, function(i) {
        s_tmp = step_matrix
        for (j in (n - m + 1):n) {
            if (i <= j) {
                s_tmp[i, j] <- (.copula_constant_joint_nll(pars + deltas[, i] + deltas[, j], spec) - g[i] - g[j] + init_llh) / s_tmp[i, j]
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

.copula_constant_scores <- function(pars, spec)
{
    estimate <- NULL
    garch_scores <- do.call(cbind, lapply(spec$univariate, function(x) .rescale_univariate_scores(x)))
    N <- NROW(garch_scores)
    n <- length(pars)
    m <- length(spec$parmatrix[estimate == 1]$value)
    .epsilon <- .Machine$double.eps
    copula_pars <- tail(pars, m)
    step_size <-  pmax(abs(copula_pars/2), 0.01) * .epsilon^(1/3)
    step_plus <- copula_pars + step_size
    step_minus <- copula_pars - step_size
    likvec_plus <- likvec_minus <- matrix(0, ncol = m, nrow = N)
    pars_plus <- pars_minus <- pars
    for (i in 1:m) {
        step_pars_plus <- copula_pars
        step_pars_minus <- copula_pars
        step_pars_plus[i] <- step_plus[i]
        step_pars_minus[i] <- step_minus[i]
        pars_plus[(n - m + 1):n] <- step_pars_plus
        pars_minus[(n - m + 1):n] <- step_pars_minus
        likvec_plus[,i] <- .copula_constant_joint_nllvec(pars_plus, spec)
        likvec_minus[,i] <- .copula_constant_joint_nllvec(pars_minus, spec)
    }
    central_diff <- likvec_plus - likvec_minus
    copula_scores <- matrix(NA, ncol = m, nrow = N)
    central_diff_sd = 2 * kronecker(matrix(1, N, 1), t(step_size))
    for (i in 1:m) {
        copula_scores[, i] = central_diff[, i]/central_diff_sd[, i]
    }
    joint_scores <- cbind(garch_scores, copula_scores)
    return(joint_scores)
}
