# !diagnostics suppress=data.table,copy,as.data.table,setnames

# copula c++ values ---------------------------------------------------

.copula_dynamic_values <- function(pars, spec, type = "nll", return_all = FALSE)
{
    estimate <- group <- NULL
    type <- match.arg(type[1], c("nll","ll_vec","Z","R","Qbar","Nbar","Q"))
    pmatrix <- copy(spec$parmatrix)
    pmatrix[estimate == 1]$value <- pars
    dist <- spec$copula
    dccorder <- spec$dynamics$order
    if (spec$dynamics$model == "adcc") {
        dccorder <- c(dccorder[1], dccorder[1], dccorder[2])
    } else {
        dccorder <- c(dccorder[1], 0, dccorder[2])
    }
    if (return_all) {
        out <- switch(dist,
                      "mvn" = .copula_dynamic_normal(alpha = pmatrix[group == "alpha"]$value,
                                                          gamma = pmatrix[group == "gamma"]$value,
                                                          beta = pmatrix[group == "beta"]$value,
                                                          u = spec$transform$u,
                                                          dccorder = dccorder),
                      "mvt" = .copula_dynamic_student(alpha = pmatrix[group == "alpha"]$value,
                                                          gamma = pmatrix[group == "gamma"]$value,
                                                          beta = pmatrix[group == "beta"]$value,
                                                          shape = pmatrix[group == "shape"]$value,
                                                          u = spec$transform$u,
                                                          dccorder = dccorder)
        )
    } else {
        out <- switch(dist,
                      "mvn" = .copula_dynamic_normal(alpha = pmatrix[group == "alpha"]$value,
                                                          gamma = pmatrix[group == "gamma"]$value,
                                                          beta = pmatrix[group == "beta"]$value,
                                                          u = spec$transform$u,
                                                          dccorder = dccorder)[[type]],
                      "mvt" = .copula_dynamic_student(alpha = pmatrix[group == "alpha"]$value,
                                                          gamma = pmatrix[group == "gamma"]$value,
                                                          beta = pmatrix[group == "beta"]$value,
                                                          shape = pmatrix[group == "shape"]$value,
                                                          u = spec$transform$u,
                                                          dccorder = dccorder)[[type]]
        )
    }
    return(out)
}

.copula_constant_values <- function(pars, spec, type = "nll", return_all = FALSE)
{
    estimate <- group <- NULL
    pmatrix <- copy(spec$parmatrix)
    if (!is.null(pars)) {
        pmatrix[estimate == 1]$value <- pars
    }
    dist <- spec$copula
    if (return_all) {
        out <- switch(dist,
                      "mvn" = .copula_constant_normal(u = spec$transform$u, method = spec$dynamics$constant),
                      "mvt" = .copula_constant_student(u = spec$transform$u, shape = pmatrix[group == "shape"]$value)
        )
    } else {
        out <- switch(dist,
                      "mvn" = .copula_constant_normal(u = spec$transform$u, method = spec$dynamics$constant)[[type]],
                      "mvt" = .copula_constant_student(u = spec$transform$u, shape = pmatrix[group == "shape"]$value)[[type]]
        )
    }
    return(out)
}

# copula c++ fun ---------------------------------------------------

.copula_dynamic_fun <- function(spec, type = "nll")
{
    estimate <- group <- NULL
    if (spec$copula == "mvn") {
        fun <- function(x, arglist)
        {
            arglist$parmatrix[estimate == 1]$value <- x
            .copula_dynamic_normal(alpha = arglist$parmatrix[group == "alpha"]$value,
                                   gamma = arglist$parmatrix[group == "gamma"]$value,
                                   beta = arglist$parmatrix[group == "beta"]$value,
                                   u = arglist$u,
                                   dccorder = arglist$dccorder)[[type]]
        }
    } else {
        fun <- function(x, arglist)
        {
            arglist$parmatrix[estimate == 1]$value <- x
            .copula_dynamic_student(alpha = arglist$parmatrix[group == "alpha"]$value,
                                   gamma = arglist$parmatrix[group == "gamma"]$value,
                                   beta = arglist$parmatrix[group == "beta"]$value,
                                   shape = arglist$parmatrix[group == "shape"]$value,
                                   u = arglist$u,
                                   dccorder = arglist$dccorder)[[type]]
        }
    }
    if (spec$dynamics$model == "adcc") {
        cons <- function(x, arglist) {
            arglist$parmatrix[estimate == 1]$value <- x
            .copula_adcc_constraint(alpha = arglist$parmatrix[group == "alpha"]$value,
                            gamma = arglist$parmatrix[group == "gamma"]$value,
                            beta = arglist$parmatrix[group == "beta"]$value,
                            shape = arglist$parmatrix[group == "shape"]$value,
                            u = arglist$u,
                            dccorder = arglist$dccorder,
                            distribution = arglist$distribution) - 0.999
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


.copula_constant_fun <- function(spec, type = "nll")
{
    estimate <- group <- NULL
    if (spec$copula == "mvt") {
        fun <- function(x, arglist)
        {
            arglist$parmatrix[estimate == 1]$value <- x
            .copula_constant_student(shape = arglist$parmatrix[group == "shape"]$value,
                                   u = arglist$u)[[type]]
        }
    } else {
        stop("\ncopula not recognized.")
    }
    grad_fun <- function(x, arglist) {
        return(grad(func = arglist$fun, x = x, arglist = arglist))
    }
    hess_fun <- function(x, arglist) {
        nderiv_hessian(func = arglist$fun, x = x, lower = arglist$lower, upper = arglist$upper, arglist = arglist)
    }
    return(list(fun = fun, grad = grad_fun, hess = hess_fun, inequality_cons = NULL, inequality_jac = NULL))
}


# copula loglik ---------------------------------------------------

.copula_dynamic_loglik <- function(pars, arglist)
{
    cgarch_env <- arglist$cgarch_env
    assign("pars", pars, envir = cgarch_env)
    if (arglist$inequality_cons(pars, arglist) >= 0) {
        nll <- as.numeric(NA)
    } else {
        nll <- arglist$fun(pars, arglist)
    }
    if (!is.finite(nll) | is.na(nll)) {
        nll <- get("nll", cgarch_env) + 0.1 * abs(get("nll", cgarch_env))
        assign("nll", nll, envir = cgarch_env)
        return(nll)
    } else {
        assign("nll", nll, envir = cgarch_env)
        return(nll)
    }
}

.copula_constant_loglik <- function(pars, arglist)
{
    cgarch_env <- arglist$cgarch_env
    nll <- arglist$fun(pars, arglist)
    if (!is.finite(nll) | is.na(nll)) {
        nll <- get("nll", cgarch_env) + 0.1 * abs(get("nll", cgarch_env))
        assign("nll", nll, envir = cgarch_env)
        return(nll)
    } else {
        assign("nll", nll, envir = cgarch_env)
        return(nll)
    }
}


# copula estimation ---------------------------------------------------

.copula_constant_estimate <- function(object, solver = "solnp", control, return_hessian = TRUE, verbose = FALSE, ...)
{
    elapsed <- Sys.time()
    group <- parameter <- estimate <- NULL
    if (object$copula == "mvn") {
        R <- .copula_constant_values(pars = NULL, spec = object, type = "R")
        hessian <- NULL
        scores <- NULL
        L <- .generate_constant_covariance(correlation = R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
        whitened_residuals <- L[["W"]]
        H <- L[["H"]]
        # reshape the correlation and covariance to save space
        H <- .lower_tri_matrix(H, diag = TRUE)
        R <- .lower_tri_matrix(R, diag = FALSE)
        pmatrix <- copy(object$parmatrix)
        garch_nll <- sum(sapply(object$univariate, function(x) -1.0 * as.numeric(logLik(x))))
        W <- .copula_constant_values(pars = NULL, spec = object, return_all = TRUE)
        copula_nll <- W[["nll"]]
        copula_residuals <- W[["Z"]]
        # quasi-likelihood
        model_nll <- garch_nll + copula_nll
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
            cgarch_env <- new.env(hash = TRUE)
            arglist$cgarch_env <- cgarch_env
            assign("nll", 1.0, envir = cgarch_env)
            arglist$parmatrix <- copy(object$parmatrix)
            arglist$u <- object$transform$u
            fnll <- .copula_constant_fun(object, type = "nll")
            arglist$fun <- fnll$fun
            arglist$grad <- fnll$grad
            arglist$hess <- fnll$hess
            arglist$distribution <- object$copula
            arglist$lower <- lower
            arglist$upper <- upper
            pars <- init_pars
            solution <- list(converged = TRUE)
            solution$conditions <- NULL
            solution$nll <- fnll$fun(init_pars, arglist)
            W <- .copula_constant_values(pars = pars, spec = object, return_all = TRUE)
            R <- W[["R"]]
            copula_residuals <- W[["Z"]]
            L <- .generate_constant_covariance(correlation = R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
            whitened_residuals <- L[["W"]]
            H <- L[["H"]]
            # reshape the correlation and covariance to save space
            H <- .lower_tri_matrix(H, diag = TRUE)
            R <- .lower_tri_matrix(R, diag = FALSE)
            pmatrix <- copy(object$parmatrix)
            garch_nll <- sum(sapply(object$univariate, function(x) -1.0 * as.numeric(logLik(x))))
            copula_nll <- solution$nll
            # quasi-likelihood
            model_nll <- garch_nll + copula_nll
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
                hessian <- .copula_constant_hessian(all_pars, new_spec)
            } else {
                hessian <- NULL
            }
            scores <- .copula_constant_scores(all_pars, new_spec)
        } else {
            estimate <- NULL
            init_pars <- object$parmatrix[group == "shape"]$value
            lower <- object$parmatrix[group == "shape"]$lower
            upper <- object$parmatrix[group == "shape"]$upper
            solver <- match.arg(solver[1], c("solnp","optim"))
            arglist <- list()
            cgarch_env <- new.env(hash = TRUE)
            arglist$cgarch_env <- cgarch_env
            assign("nll", 1.0, envir = cgarch_env)
            arglist$parmatrix <- copy(object$parmatrix)
            arglist$u <- object$transform$u
            fnll <- .copula_constant_fun(object, type = "nll")
            arglist$fun <- fnll$fun
            arglist$grad <- fnll$grad
            arglist$hess <- fnll$hess
            arglist$distribution <- object$copula
            arglist$lower <- lower
            arglist$upper <- upper
            assign("pars", init_pars, envir = cgarch_env)
            sol <- .dcc_constant_solver(solver = solver, pars = init_pars, fun = .copula_constant_loglik, lower = lower, upper = upper,
                                           control = control, arglist = arglist)
            if (inherits(sol, 'try-error')) {
                warning("\nmodel did not converge and possibly errored-out.  Returning empty list.")
                return(list())
            } else {
                pars <- .solver_extract_pars(sol, solver)
                solution <- list(converged = TRUE)
                solution$conditions <- solver_conditions(pars, .copula_constant_loglik, gr = arglist$grad, hess = arglist$hess, arglist = arglist)
                solution$nll <- .solver_extract_solution(sol, solver)
                W <- .copula_constant_values(pars = pars, spec = object, return_all = TRUE)
                R <- W[["R"]]
                copula_residuals <- W[["Z"]]
                L <- .generate_constant_covariance(correlation = R, sigmas = coredata(sigma(object$univariate)), residuals = coredata(residuals(object$univariate)))
                whitened_residuals <- L[["W"]]
                H <- L[["H"]]
                # reshape the correlation and covariance to save space
                H <- .lower_tri_matrix(H, diag = TRUE)
                R <- .lower_tri_matrix(R, diag = FALSE)
                pmatrix <- copy(object$parmatrix)
                garch_nll <- sum(sapply(object$univariate, function(x) -1.0 * as.numeric(logLik(x))))
                copula_nll <- solution$nll
                # quasi-likelihood
                model_nll <- garch_nll + copula_nll
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
                    hessian <- .copula_constant_hessian(all_pars, new_spec)
                } else {
                    hessian <- NULL
                }
                scores <- .copula_constant_scores(all_pars, new_spec)
            }
        }

    }
    elapsed <- Sys.time() - elapsed
    out <- list(parmatrix = pmatrix, solution = solution, mu = object$target$mu, hessian = hessian, scores = scores,
                R = R, H = H, whitened_residuals = whitened_residuals,
                copula_residuals = copula_residuals,
                copula_nll = copula_nll, loglik = model_nll, joint_parmatrix = mmatrix, spec = object)
    class(out) <- c("cgarch.estimate","cgarch.constant")
    return(out)
}

.copula_dynamic_estimate <- function(object, solver = "solnp", control, return_hessian = TRUE, verbose = FALSE, ...)
{
    elapsed <- Sys.time()
    estimate <- NULL
    solver <- match.arg(solver[1], c("solnp","nloptr"))
    arglist <- list()
    cgarch_env <- new.env(hash = TRUE)
    arglist$cgarch_env <- cgarch_env
    assign("nll", 1.0, envir = cgarch_env)
    arglist$parmatrix <- copy(object$parmatrix)
    arglist$u <- object$transform$u
    fnll <- .copula_dynamic_fun(object, type = "nll")
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
    arglist$distribution <- object$copula
    lower <- object$parmatrix[estimate == 1]$lower
    upper <- object$parmatrix[estimate == 1]$upper
    arglist$lower <- lower
    arglist$upper <- upper
    init_pars <- object$parmatrix[estimate == 1]$value
    maxpq <- max(object$dynamics$order)
    assign("pars", init_pars, envir = cgarch_env)
    if (verbose) {
        cat("\nestimation...", sep = "")
    }
    sol <- .dcc_dynamic_solver(solver = solver, pars = init_pars, fun = .copula_dynamic_loglik, lower = lower, upper = upper, control = control, arglist = arglist)
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
        solution$conditions <- solver_conditions(pars, .copula_dynamic_loglik, gr = arglist$grad, hess = arglist$hess, arglist = arglist)
        solution$nll <- .solver_extract_solution(sol, solver)
        W <- .copula_dynamic_values(pars = pars, spec = object, return_all = TRUE)
        R <- W[["R"]]
        Qbar <- W[["Qbar"]]
        Nbar <- W[["Nbar"]]
        Q <- W[["Q"]]
        copula_residuals <- W[["Z"]]
        if (maxpq > 0) {
            R <- R[,,-(1:maxpq)]
            copula_residuals <- copula_residuals[-(1:maxpq), , drop = FALSE]
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
        garch_nll <- sum(sapply(object$univariate, function(x) -1.0 * as.numeric(logLik(x))))
        copula_nll <- solution$nll
        # quasi-likelihood
        model_nll <- garch_nll + copula_nll
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
            hessian <- .copula_dynamic_hessian(all_pars, new_spec)
        } else {
            hessian <- NULL
        }
        scores <- .copula_dynamic_scores(all_pars, new_spec)
        elapsed <- Sys.time() - elapsed
        out <- list(parmatrix = pmatrix, solution = solution, mu = object$target$mu, hessian = hessian,
                    scores = scores, joint_parmatrix = mmatrix,
                    Qbar = Qbar, Nbar = Nbar, R = R, H = H, Q = Q,
                    copula_residuals = copula_residuals, whitened_residuals = whitened_residuals,
                    copula_nll = copula_nll, loglik = model_nll, spec = object, elapsed = elapsed)
        class(out) <- c("cgarch.estimate","cgarch.dynamic")
        return(out)
    }
}

# copula filtering ---------------------------------------------------

.copula_dynamic_filter <- function(object, y, update = TRUE, cond_mean = NULL, ...)
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
    tmp <- copula_transformation_filter(new_univariate, object$spec$transformation, object$spec$transform$transform_model, n_update = n_update)
    tmp$u[tmp$u < 3.330669e-16] <- 2.220446e-16
    tmp$u[tmp$u > 0.99999] <- 0.99999
    object$spec$transform$u <- tmp$u
    object$spec$transform$transform_model <- tmp$transform_model
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
    maxpq <- max(object$spec$dynamics$order)
    cfit <- switch(object$spec$copula,
                   "mvn" = .copula_dynamic_normal_filter(alpha = alpha, gamma = gamma, beta = beta, u = tmp$u, dccorder = dccorder, n_update = n_update),
                   "mvt" = .copula_dynamic_student_filter(alpha = alpha, gamma = gamma, beta = beta, shape = shape, u = tmp$u, dccorder = dccorder, n_update = n_update)
    )
    R <- cfit[["R"]]
    nllvec <- cfit[["ll_vec"]]
    nll <- cfit[["nll"]]
    Qbar <- cfit[["Qbar"]]
    Nbar <- cfit[["Nbar"]]
    Q <- cfit[["Q"]]
    copula_residuals <- cfit[["Z"]]
    if (maxpq > 0) {
        R <- R[,,-(1:maxpq)]
        copula_residuals <- copula_residuals[-(1:maxpq), , drop = FALSE]
        Q <- Q[,,-(1:maxpq)]
    }
    L <- .generate_dynamic_covariance(R, sigmas = coredata(sigma(new_univariate)), residuals = coredata(residuals(new_univariate)))
    whitened_residuals <- L[["W"]]
    H <- L[["H"]]
    # reshape the correlation and covariance to save space
    H <- .lower_tri_matrix(H, diag = TRUE)
    R <- .lower_tri_matrix(R, diag = FALSE)
    Q <- .lower_tri_matrix(Q, diag = TRUE)
    garch_nll <- sum(sapply(new_univariate, function(x) -1.0 * as.numeric(logLik(x))))
    # quasi-likelihood
    model_nll <- garch_nll + nll
    # hessian and scores
    object$Qbar <- Qbar
    object$Nbar <- Nbar
    object$R <- R
    object$H <- H
    object$Q <- Q
    object$copula_residuals <- copula_residuals
    object$whitened_residuals <- whitened_residuals
    object$copula_nll <- nll
    object$loglik <- model_nll
    object$elapsed <- Sys.time() - elapsed
    class(object) <- c("cgarch.estimate","cgarch.dynamic", "cgarch.filter")
    return(object)
}

.copula_constant_filter <- function(object, y, update = TRUE, cond_mean = NULL, ...)
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
    tmp <- copula_transformation_filter(new_univariate, object$spec$transformation, object$spec$transform$transform_model, n_update = n_update)
    tmp$u[tmp$u < 3.330669e-16] <- 2.220446e-16
    tmp$u[tmp$u > 0.99999] <- 0.99999
    object$spec$transform$u <- tmp$u
    object$spec$transform$transform_model <- tmp$transform_model
    shape <- object$parmatrix[group == "shape"]$value
    cfit <- switch(object$spec$copula,
                   "mvn" = .copula_constant_normal_filter(u = tmp$u, method = object$spec$dynamics$constant, n_update = n_update),
                   "mvt" = .copula_constant_student_filter(u = tmp$u, shape = shape, n_update = n_update)
    )
    R <- cfit[["R"]]
    copula_residuals <- cfit[["Z"]]
    nllvec <- cfit[["ll_vec"]]
    nll <- cfit$nll
    L <- .generate_constant_covariance(R, sigmas = coredata(sigma(new_univariate)), residuals = coredata(residuals(new_univariate)))
    whitened_residuals <- L[["W"]]
    H <- L[["H"]]
    # reshape the correlation and covariance to save space
    H <- .lower_tri_matrix(H, diag = TRUE)
    R <- .lower_tri_matrix(R, diag = FALSE)
    garch_nll <- sum(sapply(new_univariate, function(x) -1.0 * as.numeric(logLik(x))))
    # quasi-likelihood
    model_nll <- garch_nll + nll
    # hessian and scores
    object$R <- R
    object$H <- H
    object$whitened_residuals <- whitened_residuals
    object$copula_residuals <- copula_residuals
    object$copula_nll <- nll
    object$loglik <- model_nll
    object$elapsed <- Sys.time() - elapsed
    class(object) <- c("cgarch.estimate","cgarch.constant","cgarch.filter")
    return(object)
}


# copula simulation ---------------------------------------------------

.copula_dynamic_simulate_r <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0,
                                     Q_init = NULL, Z_init = NULL,
                                     init_method = c("start", "end"),
                                     cond_mean = NULL,
                                     sim_method = c("parametric", "bootstrap"), ...)
{
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    init_method <- match.arg(init_method, c("start", "end"))
    sim_method <- match.arg(sim_method, c("parametric", "bootstrap"))
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    group <- NULL
    Z <- object$copula_residuals
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
    distribution <- object$spec$copula
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
        sim <- .copula_dynamic_simulate(alpha = alpha, gamma = gamma, beta = beta, shape = shape,
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
    U_list <- lapply(sim_list, function(x) x[["U"]])
    vechn <- length(R_list[[1]][1,])
    std_residuals <- array(0, dim = c(nsim, h, n_series))
    copula_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    R <- array(unlist(R_list), dim = c(h, vechn, nsim))
    U <- array(unlist(U_list), dim = c(h, n_series, nsim))
    # for U we need the array to be concentrated on n_series not nsim
    U <- aperm(U, perm = c(3, 1, 2))
    for (i in 1:n_series) {
        std_residuals[,,i] <- .copula_qtransform(object, .retain_dimensions_array(U, i), i)
    }
    gsim <- NULL
    gsim <- future_lapply(1:n_series, function(i){
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    }, future.packages = c("tsgarch","tsmarch"), future.seed = TRUE)
    gsim <- eval(gsim)
    # reformat output [h n_series nsim]
    sim_mu <- lapply(gsim, function(x) x$series)
    sim_sigma <- lapply(gsim, function(x) x$sigma)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, n_series))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
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
    class(out) <- c("cgarch.simulate","cgarch.dynamic")
    return(out)
}

.copula_constant_simulate_r <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0, cond_mean = NULL, sim_method = "parametric", init_method = c("start", "end"), ...)
{
    # simulate from t-copula
    parameter <- NULL
    elapsed <- Sys.time()
    if (!is.null(seed)) set.seed(seed)
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    init_method <- match.arg(init_method, c("start", "end"))
    h <- h + burn
    n_series <- object$spec$n_series
    shape <- object$parmatrix[parameter == "shape"]$value
    Z <- object$copula_residuals
    R <- tscor(object)
    n_series <- object$spec$n_series
    distribution <- object$spec$copula
    sim_list <- future_lapply(1:nsim, function(i) {
        if (sim_method == "parametric") {
            std_noise <- matrix(rnorm(h * n_series), ncol = n_series, nrow = h)
        } else {
            tmp_std_noise <- .decorrelate_errors(R, Z)
            indices <- 1:NROW(tmp_std_noise)
            std_noise <- tmp_std_noise[sample(indices, h, replace = TRUE), ]
        }
        sim <- .copula_constant_simulate(shape = shape, R = R, std_noise = std_noise, timesteps = h, distribution = distribution)
        Z <- sim[["Z"]]
        U <- sim[["U"]]
        return(list(Z = Z, U = U))
    }, future.seed = TRUE, future.globals = TRUE, future.packages = "tsmarch")
    sim_list <- eval(sim_list)
    Z_list <- lapply(sim_list, function(x) x[["Z"]])
    U_list <- lapply(sim_list, function(x) x[["U"]])
    rvec <- .lower_tri_matrix(R, diag = FALSE)
    vechn <- length(rvec)
    std_residuals <- array(0, dim = c(nsim, h, n_series))
    copula_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    U <- array(unlist(U_list), dim = c(h, n_series, nsim))
    U <- aperm(U, perm = c(3, 1, 2))
    u_dims <- dim(U)
    for (i in 1:n_series) {
        std_residuals[,,i] <- .copula_qtransform(object, .retain_dimensions_array(U, i), i)
    }
    gsim <- lapply(1:n_series, function(i){
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    })
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
                model = object$spec$dynamics$model, garch_sim = gsim)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("cgarch.simulate","cgarch.constant")
    return(out)
}


# copula prediction ---------------------------------------------------

.copula_dynamic_predict <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"), forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
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
    burn <- 0
    sim_method <- match.arg(sim_method, c("parametric", "bootstrap"))
    mu <- .cond_mean_spec(cond_mean, object$spec$n_series, h, object$spec$series_names)
    group <- NULL
    Z <- object$copula_residuals
    R <- tscor(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    shape <- object$parmatrix[group == "shape"]$value
    Qbar <- object$qbar
    Nbar <- object$nbar
    dccorder <- object$spec$dynamics$order
    maxpq <- max(dccorder)
    n_series <- object$spec$n_series
    n <- NROW(object$Q)
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    distribution <- object$spec$copula
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
        sim <- .copula_dynamic_simulate(alpha = alpha, gamma = gamma, beta = beta, shape = shape,
                                        Qbar = Qbar, Nbar = Nbar,
                                        Qinit = Qinit, Zinit = Zinit,
                                        std_noise = std_noise,
                                        timesteps = h, burn = 0, dccorder = dccorder,
                                        distribution = distribution)
        return(sim)
    }, future.seed = TRUE, future.globals = TRUE, future.packages = "tsmarch")
    sim_list <- eval(sim_list)
    R_list <- lapply(sim_list, function(x) x[["R"]])
    Z_list <- lapply(sim_list, function(x) x[["Z"]])
    U_list <- lapply(sim_list, function(x) x[["U"]])
    vechn <- length(R_list[[1]][1,])
    std_residuals <- array(0, dim = c(nsim, h, n_series))
    copula_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    R <- array(unlist(R_list), dim = c(h, vechn, nsim))
    U <- array(unlist(U_list), dim = c(h, n_series, nsim))
    # for U we need the array to be concentrated on n_series not nsim
    U <- aperm(U, perm = c(3, 1, 2))
    for (i in 1:n_series) {
        std_residuals[,,i] <- .copula_qtransform(object, .retain_dimensions_array(U, i), i)
    }
    gsim <- NULL
    gsim <- lapply(1:n_series, function(i){
        .garch_simulate_model(object$spec$univariate[[i]], nsim, h, burn, .retain_dimensions_array(std_residuals, i), init_method)
    })
    # reformat output [h n_series nsim]
    sim_mu <- lapply(gsim, function(x) x$series)
    sim_sigma <- lapply(gsim, function(x) x$sigma)
    sim_mu <- array(unlist(sim_mu), dim = c(nsim, h, n_series))
    sim_mu <- aperm(sim_mu, perm = c(2, 3, 1))
    if (!is.null(cond_mean)) {
        sim_mu <- .cond_mean_inject(sim_mu, mu, recenter = TRUE)
    }
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
                model = object$spec$dynamics$model, forc_dates = forc_dates)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("cgarch.predict","cgarch.dynamic")
    return(out)
}

.copula_constant_predict <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"), forc_dates = NULL,
                                     cond_mean = NULL, seed = NULL, ...)
{
    group <- NULL
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
    Z <- object$copula_residuals
    R <- tscor(object)
    n_series <- object$spec$n_series
    distribution <- object$spec$copula
    sim_list <- future_lapply(1:nsim, function(i) {
        if (sim_method == "parametric") {
            std_noise <- matrix(rnorm(h * n_series), ncol = n_series, nrow = h)
        } else {
            tmp_std_noise <- .decorrelate_errors(R, Z)
            indices <- 1:NROW(tmp_std_noise)
            std_noise <- tmp_std_noise[sample(indices, h, replace = TRUE), ]
        }
        sim <- .copula_constant_simulate(shape = shape, R = R, std_noise = std_noise, timesteps = h, distribution = distribution)
        Z <- sim[["Z"]]
        U <- sim[["U"]]
        return(list(Z = Z, U = U))
    }, future.seed = TRUE, future.globals = TRUE, future.packages = "tsmarch")
    sim_list <- eval(sim_list)
    Z_list <- lapply(sim_list, function(x) x[["Z"]])
    U_list <- lapply(sim_list, function(x) x[["U"]])
    rvec <- .lower_tri_matrix(R, diag = FALSE)
    vechn <- length(rvec)
    std_residuals <- array(0, dim = c(nsim, h, n_series))
    copula_residuals <- array(unlist(Z_list), dim = c(h, n_series, nsim))
    U <- array(unlist(U_list), dim = c(h, n_series, nsim))
    # for U we need the array to be concentrated on n_series not nsim
    U <- aperm(U, perm = c(3, 1, 2))
    u_dims <- dim(U)
    for (i in 1:n_series) {
        std_residuals[,,i] <- .copula_qtransform(object, .retain_dimensions_array(U, i), i)
    }
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
                model = object$spec$dynamics$model, forc_dates = forc_dates)
    elapsed <- Sys.time() - elapsed
    out$elapsed <- elapsed
    class(out) <- c("cgarch.predict","cgarch.constant")
    return(out)
}

# copula newsimpact ---------------------------------------------------

.copula_newsimpact_correlation <- function(object, pair = c(1,2), ...)
{
    pair <- sort(pair)
    params <- .get_dcc_params(object)
    alpha <- params$alpha
    gamma <- params$gamma
    beta <- params$beta
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    n_series <- object$spec$n_series
    Z <- object$copula_residuals
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
    out <- list(nisurface = ni, axis = zseq, pair_names = pair_names, type = "correlation", model = "copula", dynamics = object$spec$dynamics$model)
    class(out) <- c("tsmarch.newsimpact")
    return(out)
}

.copula_newsimpact_covariance <- function(object, pair = c(1, 2), ...)
{
    pair <- sort(pair)
    params <- .get_dcc_params(object)
    alpha <- params$alpha
    gamma <- params$gamma
    beta <- params$beta
    Qbar <- object$Qbar
    Nbar <- object$Nbar
    n_series <- object$spec$n_series
    Z <- object$copula_residuals
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
    out <- list(nisurface = ni, axis = zseq, pair_names = pair_names, type = "covariance", model = "copula", dynamics = object$spec$dynamics$model)
    class(out) <- c("tsmarch.newsimpact")
    return(out)
}




