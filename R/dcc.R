# !diagnostics suppress=data.table,copy,as.data.table,setnames

# dcc loglik ---------------------------------------------------

.dcc_dynamic_loglik <- function(pars, arglist)
{
    dcc_env <- arglist$dcc_env
    assign("pars", pars, envir = dcc_env)
    if (arglist$inequality_cons(pars, arglist) >= 0) {
        nll <- as.numeric(NA)
    } else {
        nll <- arglist$fun(pars, arglist)
    }
    if (!is.finite(nll) | is.na(nll)) {
        nll <- get("nll", dcc_env) + 0.1 * abs(get("nll", dcc_env))
        assign("nll", nll, envir = dcc_env)
        return(nll)
    } else {
        assign("nll", nll, envir = dcc_env)
        return(nll)
    }
}

.dcc_constant_loglik <- function(pars, arglist)
{
    dcc_env <- arglist$dcc_env
    nll <- arglist$fun(pars, arglist)
    if (!is.finite(nll) | is.na(nll)) {
        nll <- get("nll", dcc_env) + 0.1 * abs(get("nll", dcc_env))
        assign("nll", nll, envir = dcc_env)
        return(nll)
    } else {
        assign("nll", nll, envir = dcc_env)
        return(nll)
    }
}

# dcc estimation ---------------------------------------------------

.dcc_constant_estimate <- function(object, solver = "solnp", control, return_hessian = TRUE, verbose = FALSE, ...)
{
    elapsed <- Sys.time()
    group <- parameter <- estimate <- NULL
    if (object$distribution == "mvn") {
        R <- .dcc_constant_values(pars = NULL, spec = object, type = "R")
        hessian <- NULL
        scores <- NULL
        L <- .generate_constant_covariance(correlation = R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
        whitened_residuals <- L[["W"]]
        H <- L[["H"]]
        # reshape the correlation and covariance to save space
        H <- .lower_tri_matrix(H, diag = TRUE)
        R <- .lower_tri_matrix(R, diag = FALSE)
        pmatrix <- copy(object$parmatrix)
        W <- .dcc_constant_values(pars = NULL, spec = object, return_all = TRUE)
        model_nll <- W[["nll"]]
        # quasi-likelihood
        solution <- NULL
        pmatrix <- copy(object$parmatrix)
        all_pars <- c(as.vector(sapply(object$univariate, function(x) x$parmatrix[estimate == 1]$value)))
        new_spec <- copy(object)
        new_spec$parmatrix <- copy(pmatrix)
        mmatrix <- .joint_parameter_matrix(new_spec)
    } else {
        if (object$parmatrix[group == "shape"]$estimate == 0) {
            init_pars <- object$parmatrix[group == "shape"]$value
            lower <- object$parmatrix[group == "shape"]$lower
            upper <- object$parmatrix[group == "shape"]$upper
            solver <- match.arg(solver[1], c("solnp","optim"))
            arglist <- list()
            dcc_env <- new.env(hash = TRUE)
            arglist$dcc_env <- dcc_env
            assign("nll", 1.0, envir = dcc_env)
            arglist$parmatrix <- copy(object$parmatrix)
            arglist$z <- residuals(object$univariate, standardize = TRUE)
            fnll <- .dcc_constant_fun(object, type = "nll")
            arglist$fun <- fnll$fun
            arglist$grad <- fnll$grad
            arglist$hess <- fnll$hess
            arglist$distribution <- object$distribution
            arglist$lower <- lower
            arglist$upper <- upper
            arglist$Sigma <- coredata(sigma(object$univariate))
            pars <- init_pars
            solution <- list(converged = TRUE)
            solution$conditions <- NULL
            solution$nll <- fnll$fun(init_pars, arglist)
            W <- .dcc_constant_values(pars = pars, spec = object, return_all = TRUE)
            R <- W[["R"]]
            L <- .generate_constant_covariance(correlation = R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
            whitened_residuals <- L[["W"]]
            H <- L[["H"]]
            # reshape the correlation and covariance to save space
            H <- .lower_tri_matrix(H, diag = TRUE)
            R <- .lower_tri_matrix(R, diag = FALSE)
            pmatrix <- copy(object$parmatrix)
            model_nll <- solution$nll
            # quasi-likelihood
            solution <- NULL
            pmatrix <- copy(object$parmatrix)
            # swttch to calculate
            pmatrix[parameter == "shape", estimate := 1]
            pmatrix[estimate == 1]$value <- pars
            # hessian/scores
            all_pars <- c(as.vector(sapply(object$univariate, function(x) x$parmatrix[estimate == 1]$value)), pmatrix[estimate == 1]$value)
            new_spec <- copy(object)
            new_spec$parmatrix <- copy(pmatrix)
            mmatrix <- .joint_parameter_matrix(new_spec)
            new_spec$joint_parmatrix <- mmatrix
            if (verbose) {
                cat("\nhessian and score calculation.", sep = "")
            }
            if (return_hessian) {
                hessian <- .dcc_constant_hessian(all_pars, new_spec)
            } else {
                hessian <- NULL
            }
            scores <- .dcc_constant_scores(all_pars, new_spec)
        } else {
            estimate <- NULL
            init_pars <- object$parmatrix[group == "shape"]$value
            lower <- object$parmatrix[group == "shape"]$lower
            upper <- object$parmatrix[group == "shape"]$upper
            solver <- match.arg(solver[1], c("solnp","optim"))
            arglist <- list()
            dcc_env <- new.env(hash = TRUE)
            arglist$dcc_env <- dcc_env
            assign("nll", 1.0, envir = dcc_env)
            arglist$parmatrix <- copy(object$parmatrix)
            arglist$z <- residuals(object$univariate, standardize = TRUE)
            fnll <- .dcc_constant_fun(object, type = "nll")
            arglist$fun <- fnll$fun
            arglist$grad <- fnll$grad
            arglist$hess <- fnll$hess
            arglist$distribution <- object$distribution
            arglist$lower <- lower
            arglist$upper <- upper
            arglist$s <- coredata(sigma(object$univariate))
            assign("pars", init_pars, envir = dcc_env)
            sol <- .dcc_constant_solver(solver = solver, pars = init_pars, fun = .dcc_constant_loglik, lower = lower, upper = upper,
                                           control = control, arglist = arglist)
            if (inherits(sol, 'try-error')) {
                warning("\nmodel did not converge and possibly errored-out.  Returning empty list.")
                return(list())
            } else {
                pars <- .solver_extract_pars(sol, solver)
                solution <- list(converged = TRUE)
                solution$conditions <- solver_conditions(pars, .dcc_constant_loglik, gr = arglist$grad, hess = arglist$hess, arglist = arglist)
                solution$nll <- .solver_extract_solution(sol, solver)
                W <- .dcc_constant_values(pars = pars, spec = object, return_all = TRUE)
                R <- W[["R"]]
                L <- .generate_constant_covariance(correlation = R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
                whitened_residuals <- L[["W"]]
                H <- L[["H"]]
                # reshape the correlation and covariance to save space
                H <- .lower_tri_matrix(H, diag = TRUE)
                R <- .lower_tri_matrix(R, diag = FALSE)
                pmatrix <- copy(object$parmatrix)
                model_nll <- solution$nll
                # quasi-likelihood
                solution <- NULL
                pmatrix <- copy(object$parmatrix)
                pmatrix[estimate == 1]$value <- pars
                # hessian/scores
                all_pars <- c(as.vector(sapply(object$univariate, function(x) x$parmatrix[estimate == 1]$value)), pmatrix[estimate == 1]$value)
                new_spec <- copy(object)
                new_spec$parmatrix <- copy(pmatrix)
                mmatrix <- .joint_parameter_matrix(new_spec)
                new_spec$joint_parmatrix <- mmatrix
                if (verbose) {
                    cat("\nhessian and score calculation.", sep = "")
                }
                if (return_hessian) {
                    hessian <- .dcc_constant_hessian(all_pars, new_spec)
                } else {
                    hessian <- NULL
                }
                scores <- .dcc_constant_scores(all_pars, new_spec)
            }
        }
    }
    elapsed <- Sys.time() - elapsed
    out <- list(parmatrix = pmatrix, solution = solution, mu = object$target$mu, hessian = hessian, scores = scores,
                R = R, H = H, whitened_residuals = whitened_residuals,
                loglik = model_nll, joint_parmatrix = mmatrix, spec = object)
    class(out) <- c("dcc.estimate","dcc.constant")
    return(out)
}

.dcc_dynamic_estimate <- function(object, solver = "solnp", control, return_hessian = TRUE, verbose = FALSE, ...)
{
    elapsed <- Sys.time()
    estimate <- NULL
    solver <- match.arg(solver[1], c("solnp","nloptr"))
    arglist <- list()
    dcc_env <- new.env(hash = TRUE)
    arglist$dcc_env <- dcc_env
    assign("nll", 1.0, envir = dcc_env)
    arglist$parmatrix <- copy(object$parmatrix)
    arglist$z <- residuals(object$univariate, standardize = TRUE)
    fnll <- .dcc_dynamic_fun(object, type = "nll")
    arglist$fun <- fnll$fun
    arglist$grad <- fnll$grad
    arglist$hess <- fnll$hess
    arglist$inequality_cons <- fnll$inequality_cons
    arglist$inequality_jac <- fnll$inequality_jac
    dccorder <- object$dynamics$order
    if (object$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    arglist$dccorder <- dccorder
    arglist$distribution <- object$distribution
    lower <- object$parmatrix[estimate == 1]$lower
    upper <- object$parmatrix[estimate == 1]$upper
    arglist$lower <- lower
    arglist$upper <- upper
    arglist$s <- coredata(sigma(object$univariate))
    init_pars <- object$parmatrix[estimate == 1]$value
    maxpq <- max(object$dynamics$order)
    assign("pars", init_pars, envir = dcc_env)
    if (verbose) {
        cat("\nestimation...", sep = "")
    }
    sol <- .dcc_dynamic_solver(solver = solver, pars = init_pars, fun = .dcc_dynamic_loglik, lower = lower, upper = upper, control = control, arglist = arglist)
    if (inherits(sol, 'try-error')) {
        if (verbose) {
            cat("\n...failure", sep = "")
        }
        warning("\nmodel did not converge and possibly errored-out.  Returning empty list.")
        return(list())
    } else {
        if (verbose) {
            cat("\n...success", sep = "")
        }
        pars <- .solver_extract_pars(sol, solver)
        solution <- list(converged = TRUE)
        solution$conditions <- solver_conditions(pars, .dcc_dynamic_loglik, gr = arglist$grad, hess = arglist$hess, arglist = arglist)
        solution$nll <- .solver_extract_solution(sol, solver)
        W <- .dcc_dynamic_values(pars = pars, spec = object, return_all = TRUE)
        R <- W[["R"]]
        Qbar <- W[["Qbar"]]
        Nbar <- W[["Nbar"]]
        Q <- W[["Q"]]
        if (maxpq > 0) {
            R <- R[,,-(1:maxpq)]
            Q <- Q[,,-(1:maxpq)]
        }
        L <- .generate_dynamic_covariance(R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
        whitened_residuals <- L[["W"]]
        H <- L[["H"]]
        # reshape the correlation and covariance to save space
        H <- .lower_tri_matrix(H, diag = TRUE)
        R <- .lower_tri_matrix(R, diag = FALSE)
        Q <- .lower_tri_matrix(Q, diag = TRUE)
        pmatrix <- copy(object$parmatrix)
        pmatrix[estimate == 1]$value <- pars
        dcc_nll <- W[["dcc_nll"]]
        model_nll <- solution$nll
        # hessian and scores
        all_pars <- c(as.vector(sapply(object$univariate, function(x) x$parmatrix[estimate == 1]$value)), pmatrix[estimate == 1]$value)
        new_spec <- copy(object)
        new_spec$parmatrix <- copy(pmatrix)
        mmatrix <- .joint_parameter_matrix(new_spec)
        new_spec$joint_parmatrix <- mmatrix
        if (verbose) {
            cat("\nhessian and score calculation.", sep = "")
        }
        if (return_hessian) {
            hessian <- .dcc_dynamic_hessian(all_pars, new_spec)
        } else {
            hessian <- NULL
        }
        scores <- .dcc_dynamic_scores(all_pars, new_spec)
        elapsed <- Sys.time() - elapsed
        out <- list(parmatrix = pmatrix, solution = solution, mu = object$target$mu, hessian = hessian,
                    scores = scores, joint_parmatrix = mmatrix,
                    Qbar = Qbar, Nbar = Nbar, R = R, H = H, Q = Q,
                    whitened_residuals = whitened_residuals,
                    dcc_nll = dcc_nll, loglik = model_nll, spec = object, elapsed = elapsed)
        class(out) <- c("dcc.estimate","dcc.dynamic")
        return(out)
    }
}

# dcc filtering ---------------------------------------------------

.dcc_dynamic_filter <- function(object, y, update = TRUE, cond_mean = NULL, ...)
{
    group <- NULL
    # the filtering will return a new estimate object with updated data.
    # the filtering will always use information upto to T (T<y) to update
    # the intercepts so the function can be called iteratively with 1-step
    # ahead or one shot with n-step ahead.
    # This can be switch off by using update which resets the calculation based
    # on the original data size
    elapsed <- Sys.time()
    if (!is.xts(y)) stop("\ny must be an xts object.")
    if (!is.null(y)) {
        is_null_y <- FALSE
        new_y <- NROW(y)
    } else {
        is_null_y <- TRUE
    }
    y <- .check_y_filter(object, y = y)
    # filter y
    if (update) {
        n_update <- NROW(object$spec$target$y)
    } else {
        if (!is.null(object$spec$target$original_size)) {
            n_update <- object$spec$target$original_size
        } else {
            n_update <- NROW(object$spec$target$y)
        }
    }
    if (is.null(object$spec$target$original_size)) {
        object$spec$target$original_size <- NROW(object$spec$target$y)
    }
    new_univariate <- .garch_filter_model(object, y)
    object$spec$target$y <- coredata(residuals(new_univariate))
    if (!is_null_y) {
        if (!is.null(cond_mean)) {
            mu <- .cond_mean_spec(mu = cond_mean, object$spec$n_series, new_y, object$spec$series_names)
            object$spec$target$mu <- rbind(object$mu, mu)
            object$mu <- rbind(object$mu, mu)
        } else {
            object$spec$target$mu <- coredata(fitted(new_univariate))
            object$mu <- coredata(fitted(new_univariate))
        }
    } else {
        object$spec$target$mu <- coredata(fitted(new_univariate))
        object$mu <- coredata(fitted(new_univariate))
    }
    object$spec$target$sigma <- coredata(sigma(new_univariate))
    object$spec$target$index <- .garch_extract_index(new_univariate)
    object$spec$univariate <- new_univariate
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    shape <- object$parmatrix[group == "shape"]$value
    dccorder <- object$spec$dynamics$order
    if (object$spec$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    z <- residuals(new_univariate, standardize = TRUE)
    s <- coredata(sigma(new_univariate))
    maxpq <- max(object$spec$dynamics$order)
    cfit <- switch(object$spec$distribution,
                   "mvn" = .dcc_dynamic_normal_filter(alpha = alpha, gamma = gamma, beta = beta, z = z, s = s, dccorder = dccorder, n_update = n_update),
                   "mvt" = .dcc_dynamic_student_filter(alpha = alpha, gamma = gamma, beta = beta, shape = shape, z = z, s = s, dccorder = dccorder, n_update = n_update)
    )
    R <- cfit[["R"]]
    nllvec <- cfit[["ll_vec"]]
    nll <- cfit[["nll"]]
    dcc_nll <- cfit[["dcc_nll"]]
    Qbar <- cfit[["Qbar"]]
    Nbar <- cfit[["Nbar"]]
    Q <- cfit[["Q"]]
    if (maxpq > 0) {
        R <- R[,,-(1:maxpq)]
        Q <- Q[,,-(1:maxpq)]
    }
    L <- .generate_dynamic_covariance(R, sigmas = coredata(sigma(new_univariate)), residuals = coredata(residuals(new_univariate)))
    whitened_residuals <- L[["W"]]
    H <- L[["H"]]
    # reshape the correlation and covariance to save space
    H <- .lower_tri_matrix(H, diag = TRUE)
    R <- .lower_tri_matrix(R, diag = FALSE)
    Q <- .lower_tri_matrix(Q, diag = TRUE)
    # quasi-likelihood
    model_nll <- nll
    # hessian and scores
    object$Qbar <- Qbar
    object$Nbar <- Nbar
    object$R <- R
    object$H <- H
    object$Q <- Q
    object$whitened_residuals <- whitened_residuals
    object$dcc_nll <- dcc_nll
    object$loglik <- model_nll
    object$elapsed <- Sys.time() - elapsed
    class(object) <- c("dcc.estimate","dcc.dynamic", "dcc.filter")
    return(object)
}

.dcc_constant_filter <- function(object, y, update = TRUE, cond_mean = NULL, ...)
{

    elapsed <- Sys.time()
    group <- NULL
    if (!is.xts(y)) stop("\ny must be an xts object.")
    if (!is.null(y)) {
        is_null_y <- FALSE
        new_y <- NROW(y)
    } else {
        is_null_y <- TRUE
    }
    y <- .check_y_filter(object, y = y)
    # filter y
    if (update) {
        n_update <- NROW(object$spec$target$y)
    } else {
        if (!is.null(object$spec$target$original_size)) {
            n_update <- object$spec$target$original_size
        } else {
            n_update <- NROW(object$spec$target$y)
        }
    }
    if (is.null(object$spec$target$original_size)) {
        object$spec$target$original_size <- NROW(object$spec$target$y)
    }
    new_univariate <- .garch_filter_model(object, y)
    object$spec$target$y <- coredata(residuals(new_univariate))
    if (!is_null_y) {
        if (!is.null(cond_mean)) {
            mu <- .cond_mean_spec(mu = cond_mean, object$spec$n_series, new_y, object$spec$series_names)
            object$spec$target$mu <- rbind(object$mu, mu)
            object$mu <- rbind(object$mu, mu)
        } else {
            object$spec$target$mu <- coredata(fitted(new_univariate))
            object$mu <- coredata(fitted(new_univariate))
        }
    } else {
        object$spec$target$mu <- coredata(fitted(new_univariate))
        object$mu <- coredata(fitted(new_univariate))
    }
    object$spec$target$sigma <- coredata(sigma(new_univariate))
    object$spec$target$index <- .garch_extract_index(new_univariate)
    object$spec$univariate <- new_univariate
    shape <- object$parmatrix[group == "shape"]$value
    Z <- residuals(new_univariate, standardize = TRUE)
    S <- coredata(sigma(new_univariate))
    cfit <- switch(object$spec$distribution,
                   "mvn" = .dcc_constant_normal_filter(Z = Z, S = S, n_update = n_update),
                   "mvt" = .dcc_constant_student_filter(Z = Z, S = S, shape = shape, n_update = n_update)
    )
    R <- cfit[["R"]]
    nllvec <- cfit[["ll_vec"]]
    nll <- cfit[["nll"]]
    L <- .generate_constant_covariance(R, sigmas = coredata(sigma(new_univariate)), residuals = coredata(residuals(new_univariate)))
    whitened_residuals <- L[["W"]]
    H <- L[["H"]]
    # reshape the correlation and covariance to save space
    H <- .lower_tri_matrix(H, diag = TRUE)
    R <- .lower_tri_matrix(R, diag = FALSE)
    # quasi-likelihood
    model_nll <- nll
    # hessian and scores
    object$R <- R
    object$H <- H
    object$whitened_residuals <- whitened_residuals
    object$loglik <- model_nll
    object$elapsed <- Sys.time() - elapsed
    class(object) <- c("dcc.estimate","dcc.constant","dcc.filter")
    return(object)
}


# dcc simulation ---------------------------------------------------

.dcc_dynamic_simulate_r <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0,
                                    Q_init = NULL, Z_init = NULL, init_method = c("start", "end"),
                                    cond_mean = NULL, sim_method = c("parametric", "bootstrap"), ...)
{
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    init_method <- match.arg(init_method, c("start", "end"))
    sim_method <- match.arg(sim_method, c("parametric", "bootstrap"))
    group <- NULL
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    Z <- residuals(object, standardize = TRUE)
    R <- tscor(object)
    h <- h + burn
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    shape <- object$parmatrix[group == "shape"]$value
    Qbar <- object$qbar
    Nbar <- object$nbar
    dccorder <- object$spec$dynamics$order
    maxpq <- max(dccorder)
    n_series <- object$spec$n_series
    n <- length(object$Q[1,])
    dims <- c(n_series, n_series, n)
    Q <- .tril2sym(object$Q, n_series,TRUE)
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    distribution <- object$spec$distribution
    if (object$spec$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    init_states <- .dcc_sim_dynamic_initialize_values(object, Q_init, Z_init, init_method)
    Qinit <- init_states$Qinit
    Zinit <- init_states$Zinit
    exc <- maxpq + burn
    sim_list <- NULL
    sim_list <- future_lapply(1:nsim, function(i) {
        if (sim_method == "parametric") {
            std_noise <- rbind(matrix(0, ncol = n_series, nrow = maxpq), matrix(rnorm(h * n_series), ncol = n_series, nrow = h))
        } else {
            tmp_std_noise <- .decorrelate_errors(R, Z)
            indices <- 1:NROW(tmp_std_noise)
            std_noise <- rbind(matrix(0, ncol = n_series, nrow = maxpq), tmp_std_noise[sample(indices, h, replace = TRUE), ])
        }
        sim <- .dcc_dynamic_simulate(alpha = alpha, gamma = gamma, beta = beta, shape = shape,
                                        Qbar = Qbar, Nbar = Nbar,
                                        Qinit = Qinit, Zinit = Zinit,
                                        std_noise = std_noise,
                                        timesteps = h, burn = burn, dccorder = dccorder,
                                        distribution = distribution)
        return(sim)
    }, future.seed = TRUE, future.globals = TRUE, future.packages = "tsmarch")
    sim_list <- eval(sim_list)
    R_list <- lapply(sim_list, function(x) x[["R"]])
    Z_list <- lapply(sim_list, function(x) x[["Z"]])
    vechn <- length(R_list[[1]][1,])
    std_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    std_residuals <- aperm(std_residuals, perm = c(3, 1, 2))
    R <- array(unlist(R_list), dim = c(h, vechn, nsim))
    gsim <- NULL
    gsim <- future_lapply(1:n_series, function(i){
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    }, future.packages = c("tsgarch","tsmarch"), future.seed = TRUE)
    gsim <- eval(gsim)
    # reformat output [h n_series nsim]
    sim_mu <- lapply(gsim, function(x) x$series)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, n_series))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
    sim_sigma <- lapply(gsim, function(x) x$sigma)
    sim_sigma <- array(unlist(sim_sigma), dim = c(nsim, h, n_series))
    sim_sigma <- aperm(sim_sigma, perm = c(2, 3, 1))
    m <- .lower_tri_dimension(length(R[1,,1]), diag = FALSE)
    n <- length(R[,1,1])
    dims <- c(m, m, n)
    # generate Ht
    # lower tri vector [h x lower_tri_vector x nsim]
    H <- .cor2cov(R, sim_sigma, n_series)
    # time varying distribution [mu H shape] [mu R shape]
    out <- list(mu = sim_mu, H = H, R = R, Z = std_residuals,
                n_series = n_series, nsim = nsim, h = h,
                seed = seed, series_names = names(object$spec$univariate),
                model = object$spec$dynamics$model, garch_sim = gsim)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("dcc.simulate","dcc.dynamic")
    return(out)
}

.dcc_constant_simulate_r <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0, cond_mean = NULL, init_method = c("start", "end"), ...)
{
    parameter <- NULL
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    init_method <- match.arg(init_method, c("start", "end"))
    R <- object$R
    n_series <- object$spec$n_series
    shape <- object$parmatrix[parameter == "shape"]$value
    dim <- .lower_tri_dimension(length(R), diag = FALSE)
    cdim <- c(dim, dim)
    rho <- .compose_tri_matrix(object$R, cdim, diag = FALSE)
    std_residuals <- array(0, dim = c(nsim, h, n_series))
    tmp_z <- matrix(rnorm(h * dim * nsim), ncol = dim, nrow = h * nsim)
    if (object$spec$distribution == "mvt") {
        tmp_s <-  .rmvt(rho, tmp_z, shape)
    } else {
        tmp_s <-  .rmvnorm(rho, tmp_z)
    }
    for (i in 1:n_series) {
        std_residuals[,,i] <- matrix(tmp_z[,i], ncol = h, nrow = nsim)
    }
    gsim <- NULL
    gsim <- future_lapply(1:n_series, function(i) {
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    }, future.packages = c("tsgarch","tsmarch"), future.seed = TRUE, future.conditions = character(0L))
    gsim <- eval(gsim)
    # create covariance matrix
    sim_mu <- lapply(gsim, function(x) x$series)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, n_series))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
    sim_sigma <- lapply(gsim, function(x) x$sigma)
    sim_sigma <- array(unlist(sim_sigma), dim = c(nsim, h, n_series))
    sim_sigma <- aperm(sim_sigma, perm = c(2, 3, 1))
    # generate Ht
    # lower tri vector [h x lower_tri_vector x nsim]
    Ra <- array(0, dim = c(h, length(R), nsim))
    for (i in 1:nsim) {
        Ra[,,i] <- matrix(R, ncol = length(R), nrow = h, byrow = TRUE)
    }
    H <- .cor2cov(Ra, sim_sigma, n_series)
    # time varying distribution [mu H shape] [mu R shape]
    out <- list(mu = sim_mu, H = H, R = R, Z = std_residuals,
                n_series = n_series, nsim = nsim, h = h,
                seed = seed, series_names = names(object$spec$univariate),
                model = object$spec$dynamics$model, garch_sim = gsim)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("dcc.simulate","dcc.constant")
    return(out)
}

# dcc prediction ---------------------------------------------------

.dcc_dynamic_predict <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"),
                                 forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
{
    parameter <- NULL
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    if (is.null(forc_dates)) {
        forc_dates <- .forecast_dates(forc_dates, h = h, sampling = object$spec$univariate[[1]]$spec$target$sampling,
                                      last_index = tail(object$spec$target$index,1))
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h.")
    }
    init_method <- "end"
    burn <- 0
    sim_method <- match.arg(sim_method, c("parametric", "bootstrap"))
    group <- NULL
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    Z <- residuals(object, standardize = TRUE)
    R <- tscor(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    shape <- object$parmatrix[group == "shape"]$value
    dccorder <- object$spec$dynamics$order
    maxpq <- max(dccorder)
    n_series <- object$spec$n_series
    n <- NROW(object$Q)
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    distribution <- object$spec$distribution
    if (object$spec$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    Q_init <- Z_init <- NULL
    init_states <- .dcc_sim_dynamic_initialize_values(object, Q_init, Z_init, init_method)
    Qinit <- init_states$Qinit
    Zinit <- init_states$Zinit
    exc <- maxpq
    sim_list <- NULL
    sim_list <- future_lapply(1:nsim, function(i) {
        if (sim_method == "parametric") {
            std_noise <- rbind(matrix(0, ncol = n_series, nrow = maxpq), matrix(rnorm(h * n_series), ncol = n_series, nrow = h))
        } else {
            tmp_std_noise <- .decorrelate_errors(R, Z)
            indices <- 1:NROW(tmp_std_noise)
            std_noise <- rbind(matrix(0, ncol = n_series, nrow = maxpq), tmp_std_noise[sample(indices, h, replace = TRUE), ])
        }
        sim <- .dcc_dynamic_simulate(alpha = alpha, gamma = gamma, beta = beta, shape = shape,
                                     Qbar = Qbar, Nbar = Nbar,
                                     Qinit = Qinit, Zinit = Zinit,
                                     std_noise = std_noise,
                                     timesteps = h, burn = 0, dccorder = dccorder,
                                     distribution = distribution)
        sim$std_noise <- std_noise[-c(1:maxpq), ]
        return(sim)
    }, future.seed = TRUE, future.globals = TRUE, future.packages = "tsmarch")
    sim_list <- eval(sim_list)
    R_list <- lapply(sim_list, function(x) x[["R"]])
    Z_list <- lapply(sim_list, function(x) x[["Z"]])
    std_noise <- lapply(sim_list, function(x) x[["std_noise"]])
    std_noise <-  array(unlist(std_noise), dim = c(h, n_series, nsim))
    std_noise <- aperm(std_noise, perm = c(3, 1, 2))

    vechn <- length(R_list[[1]][1,])
    std_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    std_residuals <- aperm(std_residuals, perm = c(3, 1, 2))
    R <- array(unlist(R_list), dim = c(h, vechn, nsim))
    gsim <- NULL
    gsim <- lapply(1:n_series, function(i){
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    })
    # reformat output [h n_series nsim]
    sim_mu <- lapply(gsim, function(x) x$series)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, n_series))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
    sim_sigma <- lapply(gsim, function(x) x$sigma)
    sim_sigma <- array(unlist(sim_sigma), dim = c(nsim, h, n_series))
    sim_sigma <- aperm(sim_sigma, perm = c(2, 3, 1))
    m <- .lower_tri_dimension(length(R[1,,1]), diag = FALSE)
    n <- length(R[,1,1])
    dims <- c(m, m, n)
    # generate Ht
    # lower tri vector [h x lower_tri_vector x nsim]
    H <- .cor2cov(R, sim_sigma, n_series)
    out <- list(mu = sim_mu, H = H, R = R, Z = std_residuals, std_noise = std_noise,
                n_series = n_series, nsim = nsim, h = h,
                seed = seed, series_names = names(object$spec$univariate),
                sim_sigma = sim_sigma,
                model = object$spec$dynamics$model, distribution = object$spec$distribution,
                shape = object$parmatrix[parameter == "shape"]$value,
                forc_dates = forc_dates)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("dcc.predict","dcc.dynamic")
    return(out)
}

.dcc_constant_predict <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"), forc_dates = NULL,
                                  cond_mean = NULL, seed = NULL, ...)
{
    parameter <- group <- NULL
    elapsed <- Sys.time()
    init_method <- "end"
    if (!is.null(seed)) set.seed(seed)
    if (is.null(forc_dates)) {
        forc_dates <- .forecast_dates(forc_dates, h = h, sampling = object$spec$univariate[[1]]$spec$target$sampling,
                                      last_index = tail(object$spec$target$index,1))
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h.")
    }
    shape <- object$parmatrix[group == "shape"]$value
    burn <- 0
    sim_method <- match.arg(sim_method, c("parametric", "bootstrap"))
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    group <- NULL
    Z <- residuals(object, standardize = TRUE)
    R <- tscor(object)
    n_series <- object$spec$n_series
    distribution <- object$spec$distribution
    sim_list <- future_lapply(1:nsim, function(i) {
        if (sim_method == "parametric") {
            std_noise <- matrix(rnorm(h * n_series), ncol = n_series, nrow = h)
        } else {
            tmp_std_noise <- .decorrelate_errors(R, Z)
            indices <- 1:NROW(tmp_std_noise)
            std_noise <- tmp_std_noise[sample(indices, h, replace = TRUE), ]
        }
        sim <- .dcc_constant_simulate(shape = shape, R = R, std_noise = std_noise, timesteps = h, distribution = distribution)
        Z <- sim[["Z"]]
        return(list(Z = Z))
    }, future.seed = TRUE, future.globals = TRUE, future.packages = "tsmarch")
    sim_list <- eval(sim_list)
    Z_list <- lapply(sim_list, function(x) x[["Z"]])
    rvec <- .lower_tri_matrix(R, diag = FALSE)
    vechn <- length(rvec)
    std_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    std_residuals <- aperm(std_residuals, perm = c(3, 1, 2))
    # for U we need the array to be concentrated on n_series not nsim
    gsim <- lapply(1:n_series, function(i){
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    })
    # reformat output [h n_series nsim]
    sim_mu <- lapply(gsim, function(x) x$series)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, n_series))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
    sim_sigma <- lapply(gsim, function(x) x$sigma)
    sim_sigma <- array(unlist(sim_sigma), dim = c(nsim, h, n_series))
    sim_sigma <- aperm(sim_sigma, perm = c(2, 3, 1))
    # generate Ht
    # lower tri vector [h x lower_tri_vector x nsim]
    H <- .cor2cov2(matrix(rvec, nrow = 1), sim_sigma, n_series)
    # time varying distribution [mu H shape] [mu R shape]
    R <- .lower_tri_matrix(R, diag = FALSE)
    out <- list(mu = sim_mu, H = H, R = R, Z = std_residuals,
                n_series = n_series, nsim = nsim, h = h,
                seed = seed, series_names = names(object$spec$univariate),
                model = object$spec$dynamics$model, distribution = object$spec$distribution,
                shape = object$parmatrix[parameter == "shape"]$value,
                forc_dates = forc_dates)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("dcc.predict","dcc.constant")
    return(out)
}


# dcc newsimpact ---------------------------------------------------

.dcc_newsimpact_correlation <- function(object, pair = c(1,2), ...)
{
    pair <- sort(pair)
    params <- .get_dcc_params(object)
    alpha <- params$alpha
    gamma <- params$gamma
    beta <- params$beta
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    n_series <- object$spec$n_series
    Z <- residuals(object, standardize = TRUE)
    pair_names <- c(object$spec$series_names[pair[1]], object$spec$series_names[pair[2]])
    maxz <- round(max(apply(Z, 1, "max")) + 1, 0)
    minz <- round(min(apply(Z, 1, "min")) - 1, 0)
    zseq <- seq(minz, maxz, length.out = 100)
    U <- Qbar * (1 - sum(alpha) - sum(beta)) - gamma * Nbar
    ni <- matrix(0, 100, 100)
    for (i in 1:100) {
        for (j in 1:100) {
            z <- za <- matrix(0, ncol = n_series, nrow = 1)
            z[1,pair[1]] <- zseq[i]
            z[1,pair[2]] <- zseq[j]
            assymetric_z <- z * .sign_matrix(z)
            tmp <- U + alpha * t(z) %*% z + gamma * t(assymetric_z) %*% assymetric_z + beta * Qbar
            tmp <- tmp/(sqrt(diag(tmp)) %*% t(sqrt(diag(tmp))))
            ni[i,j] = tmp[pair[1],pair[2]]
        }
    }
    out <- list(nisurface = ni, axis = zseq, pair_names = pair_names, type = "correlation", model = "dcc", dynamics = object$spec$dynamics$model)
    class(out) <- c("tsmarch.newsimpact")
    return(out)
}

.dcc_newsimpact_covariance <- function(object, pair = c(1, 2), ...)
{
    pair <- sort(pair)
    params <- .get_dcc_params(object)
    alpha <- params$alpha
    gamma <- params$gamma
    beta <- params$beta
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    n_series <- object$spec$n_series
    Z <- residuals(object, standardize = TRUE)
    pair_names <- c(object$spec$series_names[pair[1]], object$spec$series_names[pair[2]])
    maxz <- round(max(apply(Z, 1, "max")) + 1, 0)
    minz <- round(min(apply(Z, 1, "min")) - 1, 0)
    zseq <- seq(minz, maxz, length.out = 100)
    unc_sigma <- sqrt(c(unconditional(object$spec$univariate[[pair[1]]]), unconditional(object$spec$univariate[[pair[2]]])))
    U <- Qbar * (1 - sum(alpha) - sum(beta)) - gamma * Nbar
    ni <- matrix(0, 100, 100)
    for (i in 1:100) {
        for (j in 1:100) {
            z <- za <- matrix(0, ncol = n_series, nrow = 1)
            z[1,pair[1]] <- zseq[i]
            z[1,pair[2]] <- zseq[j]
            assymetric_z <- z * .sign_matrix(z)
            tmp <- U + alpha * t(z) %*% z + gamma * t(assymetric_z) %*% assymetric_z + beta * Qbar
            tmp <- tmp/(sqrt(diag(tmp)) %*% t(sqrt(diag(tmp))))
            rtmp <- tmp[pair[1],pair[2]]
            ni[i,j] <- prod(unc_sigma) * rtmp
        }
    }
    out <- list(nisurface = ni, axis = zseq, pair_names = pair_names, type = "covariance", model = "dcc", dynamics = object$spec$dynamics$model)
    class(out) <- c("tsmarch.newsimpact")
    return(out)
}


.dcc_news_impact_surface <- function(x, ...)
{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    axis_names <- c(paste0(x$pair_names[1]," [t-1]"), paste0(x$pair_names[2]," [t-1]"))
    x1 <- drapecol(x$nisurface, col = topo.colors(100), NAcol = "white")
    par(mar = c(1.5,1.5,2,1.5))
    persp(  x = x$axis,
            y = x$axis,
            z = x$nisurface,  col = x1, theta = 45, phi = 25, expand = 0.5,
            ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = axis_names[1],
            ylab = axis_names[2], zlab = x$type,
            cex.axis = 0.7, cex.main = 0.8, main = paste0(toupper(x$model)," [",toupper(x$dynamics),"] News Impact ", upper_case_i(x$type,1,1)," Surface"))
    return(invisible(x))
}

# dcc predict analytic ---------------------------------------------------

.dcc_dynamic_predict_analytic <- function(object, h = 1, forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
{
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    if (is.null(forc_dates)) {
        forc_dates <- .forecast_dates(forc_dates, h = h, sampling = object$spec$univariate[[1]]$spec$target$sampling,
                                      last_index = tail(object$spec$target$index,1))
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h.")
    }
    init_method <- "end"
    group <- NULL
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    Z <- residuals(object, standardize = TRUE)
    R <- tscor(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    shape <- object$parmatrix[group == "shape"]$value
    dcc_order <- object$spec$dynamics$order
    maxpq <- max(dcc_order)
    n_series <- object$spec$n_series
    n <- NROW(object$Q)
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    distribution <- object$spec$distribution
    if (object$spec$dynamics$model == "adcc") {
        dcc_order <- c(dcc_order[1], dcc_order[1], dcc_order[2])
    } else {
        dcc_order <- c(dcc_order[1], 0, dcc_order[2])
    }
    Q_init <- Z_init <- NULL
    init_states <- .dcc_sim_dynamic_initialize_values(object, Q_init, Z_init, init_method)
    Q_init <- init_states$Qinit
    Z_init <- init_states$Zinit
    R_init <- tail(object$R, maxpq)
    H_init <- tail(object$H, maxpq)
    exc <- maxpq

    garch_predictions <- lapply(1:n_series, function(i){
        predict(object$spec$univariate[[i]], h = h, nsim  = 1)
    })

    sigmas <- do.call(cbind, lapply(1:n_series, function(i) {
        coredata(garch_predictions[[i]]$sigma)
    }))
    sigmas <- rbind(matrix(0, ncol = n_series, nrow = maxpq), sigmas)
    Z <- coredata(residuals(object, standardize = TRUE))
    s_matrix <- .sign_matrix(Z)
    if (dcc_order[2] > 0) {
        asy_Z <- s_matrix * Z
    } else {
        asy_Z <- s_matrix * 0
    }
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    sum_alpha <- sum(alpha)
    sum_beta <- sum(beta)
    sum_gamma <- sum(gamma)
    dcc_sum <- sum_alpha + sum_beta
    Omega <- Qbar * (1.0 - dcc_sum)
    if (dcc_order[1] > 0) {
        Omega <- Omega - sum_gamma * Nbar
    }
    R_predict <- H_predict <- Q_predict <- array(NA, dim = c(n_series, n_series, h + maxpq))
    for (i in 1:maxpq) {
        Q_predict[,,i] <- Q_init[,,i]
        R_predict[,,i] <- .tril2sym(matrix(R_init[i,], nrow = 1), n_series, FALSE)[,,1]
        H_predict[,,i] <- .tril2sym(matrix(H_init[i,], nrow = 1), n_series, TRUE)[,,1]
    }

    for (i in (maxpq + 1):(h + maxpq)) {
        Q_predict[,,i] <- Omega
        if (i == (maxpq + 1)) {
            if (dcc_order[1] > 0) {
                for (j in 1:dcc_order[1]) {
                    Q_predict[,,i] <- Q_predict[,,i] + alpha[j] * (Z[i - j, ] %*% t(Z[i - j, ]))
                }
            }
            if (dcc_order[2] > 0) {
                for (j in 1:dcc_order[2]) {
                    Q_predict[,,i] <- Q_predict[,,i] + gamma[j] * (asy_Z[i - j, ] %*% t(asy_Z[i - j, ]))
                }
            }
            if (dcc_order[3] > 0) {
                for (j in 1:dcc_order[3]) {
                    Q_predict[,,i] <- Q_predict[,,i] + beta[j] * Q_predict[,,i - j]
                }
            }
            Q_tmp <- diag(1/sqrt(diag(Q_predict[,,i])), n_series, n_series)
            R_predict[,,i] <- Q_tmp %*% Q_predict[,,i] %*% t(Q_tmp)
            D_tmp <- diag(sigmas[i, ], n_series, n_series)
            H_predict[,,i] <- D_tmp %*% R_predict[,,i] %*% D_tmp
            # now the unconditional calculations
            Q_star <- diag(1/sqrt(diag(Q_predict[,,i])), n_series, n_series)
            e_R <- Q_star %*% Q_predict[,,i] %*% t(Q_star)
            Q_bar_star <- diag(1/sqrt(diag(Qbar)), n_series, n_series)
            R_bar <- Q_bar_star %*% Qbar %*% t(Q_bar_star)
        } else {
            R_predict[,,i] = (1 - dcc_sum^(i - maxpq - 1)) * R_bar + dcc_sum^(i - maxpq - 1) * e_R
            D_tmp <- diag(sigmas[i, ], n_series, n_series)
            H_predict[,,i] <- D_tmp %*% R_predict[,,i] %*% D_tmp
            Q_predict[,,i] <- Q_predict[,,maxpq + 1]
        }
    }
    if (maxpq > 0) {
        Q_predict <- Q_predict[,,(maxpq + 1):(maxpq + h)]
        H_predict <- H_predict[,,(maxpq + 1):(maxpq + h)]
        R_predict <- R_predict[,,(maxpq + 1):(maxpq + h)]
    }
    out <- list(mu = cond_mean, H = H_predict, R = R_predict, Z = Z,
                n_series = n_series, h = h,
                seed = seed, series_names = names(object$spec$univariate),
                model = object$spec$dynamics$model, forc_dates = forc_dates)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("dcc.predict","dcc.dynamic")
    return(out)
}
