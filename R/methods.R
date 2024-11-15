# model specification ---------------------------------------------------


#' Generic Methods for the Copula GARCH model specification
#'
#' @param object a valid object.
#' @param ... additional parameters specific to the method implemented.
#' @returns an object of some class.
#' @aliases cgarch_modelspec
#' @rdname cgarch_modelspec.generic
#' @export
#'
cgarch_modelspec <- function(object, ...)
{
    UseMethod('cgarch_modelspec')
}

#' Copula GARCH model specification
#'
#' @param object an object of class \dQuote{tsgarch.multi_estimate}.
#' @param dynamics the type of correlation dynamics to use.
#' @param dcc_order the order of the dynamics in case of \dQuote{dcc} or \dQuote{adcc}
#' correlation models.
#' @param copula the copula distribution.
#' @param transformation the copula transformation to apply.
#' @param constant_correlation the constant correlation estimator to use. In the
#' case of the \dQuote{mvt} copula, only Kendall's tau is a valid choice.
#' @param cond_mean an optional matrix of the conditional mean for the series.
#' @param ... additional arguments passed to the \code{\link[tsdistributions]{spd_modelspec}}
#' function in the case of the \dQuote{spd} transformation.
#' @returns an object of class \dQuote{cgarch.spec}.
#' @method cgarch_modelspec tsgarch.multi_estimate
#' @rdname cgarch_modelspec
#' @export
#'
cgarch_modelspec.tsgarch.multi_estimate <- function(object, dynamics = c("constant","dcc","adcc"),
                                        dcc_order = c(1,1),
                                        transformation = c("parametric","empirical","spd"),
                                        copula = c("mvn","mvt"),
                                        constant_correlation = c("pearson", "kendall","spearman"),
                                        cond_mean = NULL, ...)
{
    dynamics <- match.arg(dynamics[1], c("constant","dcc","adcc"))
    transformation <- match.arg(transformation[1], c("parametric","empirical","spd"))
    copula <- match.arg(copula[1], c("mvn","mvt"))
    spec <- list()
    spec$target$y <- coredata(residuals(object))
    if (!is.null(cond_mean)) {
        mu <- .cond_mean_spec(mu = cond_mean, NCOL(spec$target$y), NROW(spec$target$y), names(object))
        spec$target$mu <- mu
    } else {
        spec$target$mu <- coredata(fitted(object))
    }
    spec$target$sigma <- coredata(sigma(object))
    spec$target$index <- .garch_extract_index(object)
    spec$univariate <- object
    spec$dynamics <- list()
    if (dynamics == "constant") {
        dcc_order <- c(1,1)
    } else {
        dcc_order <- c(max(0, dcc_order[1]), max(0, dcc_order[2]))
    }
    spec$dynamics$model <- dynamics
    spec$dynamics$order <- dcc_order
    spec$dynamics$constant <- match.arg(constant_correlation[1], c("pearson", "kendall","spearman"))
    if (copula == "mvt") {
        if (spec$dynamics$constant != "kendall") {
            spec$dynamics$constant = "kendall"
            warnings("\nstudent copula only available with kendall's tau correlation. Forcing to correlation kendall.")
        }
    }
    spec$transformation <- transformation
    spec$copula <- copula
    spec$parmatrix <- .copula_parameters(dynamics, copula, dcc_order)
    # add transformation here
    spec$transform <- list()
    tmp <- copula_transformation_estimate(object, transformation)
    spec$transform$u <- tmp$u
    spec$transform$transform_model <- tmp$transform_model
    # This is important else will lead to problems
    # this is what rmgarch uses:
    spec$transform$u[spec$transform$u < 3.330669e-16] <- 2.220446e-16
    spec$transform$u[spec$transform$u > 0.99999] <- 0.99999
    # number of parameters
    garch_n <- sum(sapply(spec$univariate, function(x) x$npars))
    copula_n <- sum(spec$parmatrix$estimate)
    transformation_n <- 0
    if (spec$transformation == "spd") {
        transformation_n <- 4 * length(spec$univariate)
    }
    npars <- garch_n + copula_n + transformation_n
    nobs <- NROW(spec$transform$u)
    n_series <- NCOL(spec$transform$u)
    spec$series_names <- names(object)
    spec$npars <- npars
    spec$nobs <- nobs
    spec$n_series <- n_series
    class(spec) <- "cgarch.spec"
    return(spec)
}



#' Generic Methods for the DCC GARCH model specification
#'
#' @param object a valid object.
#' @param ... additional parameters specific to the method implemented.
#' @returns an object of some class.
#' @aliases dcc_modelspec
#' @rdname dcc_modelspec.generic
#' @export
#'
dcc_modelspec <- function(object, ...)
{
    UseMethod('dcc_modelspec')
}

#' DCC GARCH model specification
#'
#' @param object an object of class \dQuote{tsgarch.multi_estimate}.
#' @param dynamics the type of correlation dynamics to use.
#' @param dcc_order the order of the dynamics in case of \dQuote{dcc} or \dQuote{adcc}
#' correlation models.
#' @param distribution the multivariate distribution. If using the \dQuote{mvt},
#' then the first stage univariate models should use the normal distribution.
#' @param cond_mean an optional matrix of the conditional mean for the series.
#' @param ... not currently used.
#' @returns an object of class \dQuote{dcc.spec}.
#' @method dcc_modelspec tsgarch.multi_estimate
#' @rdname dcc_modelspec
#' @export
#'
dcc_modelspec.tsgarch.multi_estimate <- function(object, dynamics = c("constant","dcc","adcc"),
                                        dcc_order = c(1,1),
                                        distribution = c("mvn","mvt"),
                                        cond_mean = NULL, ...)
{
    dynamics <- match.arg(dynamics[1], c("constant","dcc","adcc"))
    distribution <- match.arg(distribution[1], c("mvn","mvt"))
    spec <- list()
    spec$target$y <- coredata(residuals(object))
    if (!is.null(cond_mean)) {
        mu <- .cond_mean_spec(mu = cond_mean, NCOL(spec$target$y), NROW(spec$target$y), names(object))
        spec$target$mu <- mu
    } else {
        spec$target$mu <- coredata(fitted(object))
    }
    spec$target$sigma <- coredata(sigma(object))
    spec$target$index <- .garch_extract_index(object)
    spec$univariate <- object
    spec$dynamics <- list()
    if (dynamics == "constant") {
        dcc_order <- c(1,1)
    } else {
        dcc_order <- c(max(0, dcc_order[1]), max(0, dcc_order[2]))
    }
    spec$dynamics$model <- dynamics
    spec$dynamics$order <- dcc_order
    spec$distribution <- distribution
    spec$parmatrix <- .copula_parameters(dynamics, distribution, dcc_order)
    garch_n <- sum(sapply(spec$univariate, function(x) x$npars))
    dcc_n <- sum(spec$parmatrix$estimate)
    npars <- garch_n + dcc_n
    nobs <- NROW(spec$target$y)
    n_series <- NCOL(spec$target$y)
    spec$series_names <- names(object)
    spec$npars <- npars
    spec$nobs <- nobs
    spec$n_series <- n_series
    class(spec) <- "dcc.spec"
    return(spec)
}

# estimate ---------------------------------------------------

#' Estimates a model given a specification.
#' @description
#' Method for estimating one of the models in the package given a specification object.
#' @param object an object of class \dQuote{cgarch.spec}, \dQuote{dcc.spec} or
#' \dQuote{gogarch.spec}.
#' @param solver the solver to use for the second stage estimation of the multivariate
#' dynamics. Valid choices are \code{\link[Rsolnp]{solnp}}, \code{\link[nloptr]{nloptr}}
#' and \code{\link[stats]{optim}}, with the latter using the \dQuote{L-BFGS-B} method.
#' This option is inactive for the GOGARCH model which uses the default solver in
#' package \dQuote{tsgarch} for the estimation of the independent factors.
#' @param control solver control parameters.
#' @param return_hessian whether to calculate and return the partitioned hessian
#' of the model (see details). Not applicable in the case of the GOGARCH model.
#' @param trace whether to print tracing information for the GOGARCH model estimation.
#' @param ... for the GOGARCH model, additional options passed to the \code{\link{radical}}
#' function.
#' @returns An estimated object of one of either \dQuote{cgarch.estimate},
#' \dQuote{dcc.estimate} or \dQuote{gogarch.estimate}.
#' @details
#' ## DCC and Copula Models
#' Estimation assumes a 2-stage approach whereby the pre-estimated
#' GARCH models (first stage) are used to estimate the Copula or DCC dynamics. In
#' the case of the constant correlation Gaussian model, there are no parameters to estimate.
#' Whilst this 2-stage approach results in a very fast estimation, the calculation
#' of the standard errors based on the partitioned hessian is quite expensive. The
#' user can create a futures \code{\link[future]{plan}} to take advantage of parallel
#' functionality which for large dimensional problems leads to a large speedup.
#' Additionally, the option of not calculating the hessian (\dQuote{return_hessian})
#' is available. In that case, the scores are still calculated and the resulting
#' standard errors will be based on the outer product of gradients.
#' ## GOGARCH Model
#' The independent factors are calculated by first pre-whitening the data (PCA),
#' selecting the number of factors, and then solving for the rotation matrix. A
#' GARCH model is then estimated on each factor. A minimal amount of information
#' is retained in the estimation object and most of the heavy lifting related to
#' co-moment matrices, weighted moments etc is done through dedicated methods.
#' The estimation method makes use of parallel processing for the independent
#' factor GARCH models which can be setup using \code{\link[future]{plan}}.
#' Additionally, the RADICAL algorithm benefits from part parallelization which
#' can be controlled using \code{\link[RcppParallel]{setThreadOptions}}.
#' @author Alexios Galanos
#' @aliases estimate
#' @method estimate cgarch.spec
#' @rdname estimate.tsmarch
#' @export
#'
estimate.cgarch.spec <- function(object, solver = "solnp", control = list(trace = 0), return_hessian = TRUE, ...)
{
    control <- .default_options(control = control, solver = solver)
    out <- switch(object$dynamics$model,
           "constant" = .copula_constant_estimate(object, solver = solver, control = control, return_hessian = return_hessian, ...),
           "dcc" = .copula_dynamic_estimate(object, solver = solver, control = control, return_hessian = return_hessian, ...),
           "adcc" = .copula_dynamic_estimate(object, solver = solver, control = control, return_hessian = return_hessian, ...))
    return(out)
}


#' @method estimate dcc.spec
#' @rdname estimate.tsmarch
#' @export
#'
estimate.dcc.spec <- function(object, solver = "solnp", control = list(trace = 0), return_hessian = TRUE, ...)
{
    control <- .default_options(control = control, solver = solver)
    out <- switch(object$dynamics$model,
                  "constant" = .dcc_constant_estimate(object, solver = solver, control = control, return_hessian = return_hessian, ...),
                  "dcc" = .dcc_dynamic_estimate(object, solver = solver, control = control, return_hessian = return_hessian, ...),
                  "adcc" = .dcc_dynamic_estimate(object, solver = solver, control = control, return_hessian = return_hessian, ...))
    return(out)
}

#' @method estimate gogarch.spec
#' @rdname estimate.tsmarch
#' @export
#'
estimate.gogarch.spec <- function(object, trace = FALSE, ...)
{
    .gogarch_estimate(object, trace, ...)
}

# tsfilter ---------------------------------------------------


#' Model Filtering
#'
#' @description Filters new data based on an already estimated model.
#' @param object an object of class \dQuote{cgarch.estimate} or \dQuote{dcc.estimate}.
#' @param y an xts matrix of new values to filter.
#' @param newxreg not used in these models.
#' @param cond_mean an optional matrix of the filtered conditional mean values.
#' @param update whether to update certain values using the most recent information
#' less than the new data (see details).
#' @param ... additional arguments for future expansion.
#' @returns A \dQuote{cgarch.estimate} or \dQuote{dcc.estimate} object with updated
#' information. All values in the object are updated with the exception of the hessian
#' and scores which remain at their estimation set values.
#' @details The method filters new data and updates the object with this new information
#' so that it can be called recursively as new data arrives. The \dQuote{update} argument
#' allows values such as the intercept matrices and transformation estimates (for the
#' \dQuote{spd} and \dQuote{empirical} methods) in the dynamic case, and the constant
#' correlation in the constant case, to use information up to and include time T, where
#' T is the time stamp just preceding the new y timestamps. In this way, the filter
#' method can be called recursively and the user can selectively choose to either use the
#' updating scheme or use the original estimated values. Whatever the case, this ensures
#' that there is no look-ahead bias when filtering new information.
#' @aliases tsfilter
#' @method tsfilter cgarch.estimate
#' @rdname tsfilter.tsmarch
#' @export
#'
#'
tsfilter.cgarch.estimate <- function(object, y = NULL, newxreg = NULL, update = TRUE, cond_mean = NULL, ...)
{
    out <- switch(object$spec$dynamics$model,
                  "constant" = .copula_constant_filter(object, y = y, update = update, cond_mean = cond_mean, ...),
                  "dcc" = .copula_dynamic_filter(object, y = y, update = update, cond_mean = cond_mean, ...),
                  "adcc" = .copula_dynamic_filter(object, y = y, update = update, cond_mean = cond_mean, ...))
    return(out)
}


#' @method tsfilter dcc.estimate
#' @rdname tsfilter.tsmarch
#' @export
#'
#'
tsfilter.dcc.estimate <- function(object, y = NULL, newxreg = NULL, update = TRUE, cond_mean = NULL, ...)
{
    out <- switch(object$spec$dynamics$model,
                  "constant" = .dcc_constant_filter(object, y = y, update = update, cond_mean = cond_mean, ...),
                  "dcc" = .dcc_dynamic_filter(object, y = y, update = update, cond_mean = cond_mean, ...),
                  "adcc" = .dcc_dynamic_filter(object, y = y, update = update, cond_mean = cond_mean, ...))
    return(out)
}


#' @method tsfilter gogarch.estimate
#' @rdname tsfilter.tsmarch
#' @export
#'
tsfilter.gogarch.estimate <- function(object, y = NULL,  newxreg = NULL, cond_mean = NULL, ...)
{
    .gogarch_filter(object, y = y, cond_mean = cond_mean, ...)
}


# simulate ---------------------------------------------------

#' Model Simulation
#'
#' @description Simulates paths of a model.
#' @param object an estimated object from one of the models in the package.
#' @param nsim the number of sample paths to generate.
#' @param seed an integer that will be used in a call to set.seed before simulating.
#' @param h the number of time steps to simulate paths for.
#' @param Q_init an optional array of matrices of dimension
#' n_series x n_series x maxpq for initializing the DCC model (not relevant in
#' the constant correlation case), where maxpq is the maximum DCC model order.
#' @param Z_init an optional matrix of size maxpq x m of initialization values for the
#' standardized innovations of the DCC model. For this copula model, care should be
#' taken as these represent the DCC copula standardized innovations, not the
#' univariate GARCH innovations.
#' @param burn burn in. Will be discarded before returning the output.
#' @param init_method method to initialize the DCC and GARCH recursion (unless
#' \dQuote{Q_init} and \dQuote{Z_init} are not NULL in which case those take
#' priority for those inputs). The \dQuote{start} method initializes the recursion
#' with the same values used during estimation, whereas the \dQuote{end} method
#' uses the last values of the estimated model to initialize the recursion. In
#' the constant correlation case, only the the GARCH initialization is relevant.
#' @param sim_method white noise method for generating random sample for the
#' multivariate distribution. The default \dQuote{parametric} samples random
#' variates from the underlying error distribution whilst the \dQuote{bootstrap}
#' samples from the whitened innovations of the fitted model.
#' @param cond_mean an optional matrix (h x n_series) of the simulated conditional
#' mean for the series which is used to recenter the simulated distribution.
#' @param ... no additional arguments currently supported.
#' @details
#' Part of the code makes use of parallel functionality via
#' the \dQuote{future} package (see \code{\link[future]{plan}}). The dimension
#' the parallel execution operates on is the number of series (for the
#' individual GARCH series simulation), so unless you have more than 100 series
#' then it is possible that using a parallel back-end may actually result in
#' slower execution due to the overhead involved.
#' @returns A simulation class object for which methods exists for extracting
#' relevant statistics such as the correlation, covariance, etc.
#' @aliases simulate
#' @method simulate cgarch.estimate
#' @rdname simulate.tsmarch
#' @export
#'
#'
simulate.cgarch.estimate <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0,
                                     Q_init = NULL, Z_init = NULL,
                                     init_method = c("start", "end"), cond_mean = NULL,
                                     sim_method = c("parametric", "bootstrap"), ...)
{
    out <- switch(object$spec$dynamics$model,
                  "constant" = .copula_constant_simulate_r(object, nsim = nsim, seed = seed, h = h, burn = burn, init_method = init_method[1], sim_method = sim_method[1], cond_mean = cond_mean, ...),
                  "dcc" = .copula_dynamic_simulate_r(object, nsim = nsim, seed = seed, h = h, burn = burn, init_method = init_method[1], cond_mean = cond_mean, sim_method = sim_method[1], ...),
                  "adcc" = .copula_dynamic_simulate_r(object, nsim = nsim, seed = seed, h = h, burn = burn, init_method = init_method[1], cond_mean = cond_mean, sim_method = sim_method[1], ...))
    return(out)
}

#' @method simulate dcc.estimate
#' @rdname simulate.tsmarch
#' @export
#'
#'
simulate.dcc.estimate <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0,
                                  Q_init = NULL, Z_init = NULL,
                                  init_method = c("start", "end"), cond_mean = NULL,
                                  sim_method = c("parametric", "bootstrap"), ...)
{
    out <- switch(object$spec$dynamics$model,
                  "constant" = .dcc_constant_simulate_r(object, nsim = nsim, seed = seed, h = h, burn = burn, init_method = init_method[1], sim_method = sim_method[1], cond_mean = cond_mean, ...),
                  "dcc" = .dcc_dynamic_simulate_r(object, nsim = nsim, seed = seed, h = h, burn = burn, init_method = init_method[1], cond_mean = cond_mean, sim_method = sim_method[1], ...),
                  "adcc" = .dcc_dynamic_simulate_r(object, nsim = nsim, seed = seed, h = h, burn = burn, init_method = init_method[1], cond_mean = cond_mean, sim_method = sim_method[1], ...))
    return(out)
}

#' @method simulate gogarch.estimate
#' @rdname simulate.tsmarch
#' @export
#'
#'
simulate.gogarch.estimate <- function(object, nsim = 1, seed = NULL, h = 100, burn = 0, cond_mean = NULL, sim_method = c("parametric", "bootstrap"), ...)
{
    out <- .gogarch_simulate(object, nsim = nsim, seed = seed, h = h, burn = burn, cond_mean = cond_mean, sim_method = sim_method, ...)
    return(out)
}


# predict ---------------------------------------------------

#' Model Prediction
#'
#' @description Prediction function for estimated objects.
#' @param object an estimated object from one of the models in the package.
#' @param h the forecast horizon.
#' @param nsim the number of sample paths to generate.
#' @param seed an integer that will be used in a call to set.seed before simulating.
#' @param sim_method white noise method for generating random sample for the
#' multivariate distribution. The default \dQuote{parametric} samples random normal
#' variates whilst the \dQuote{bootstrap} samples from the whitened innovations
#' of the fitted model.
#' @param forc_dates an optional vector of forecast dates equal to h. If NULL will
#' use the implied periodicity of the data to generate a regular sequence of dates
#' after the last available date in the data.
#' @param cond_mean an optional matrix (h x n_series) of the predicted conditional
#' mean for the series which is used to recenter the simulated predictive distribution.
#' @param ... no additional arguments currently supported.
#' @details
#' For the Copula GARCH model, the prediction is based on simulation due to the
#' nonlinear transformation present in the model.
#' @returns A prediction class object for which methods exists for extracting
#' relevant statistics such as the correlation, covariance, etc.
#' @aliases predict
#' @method predict cgarch.estimate
#' @rdname predict.tsmarch
#' @export
#'
#'
predict.cgarch.estimate <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"),
                                    forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
{
    out <- switch(object$spec$dynamics$model,
                  "constant" = .copula_constant_predict(object, h = h, nsim = nsim, sim_method = sim_method[1], forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...),
                  "dcc" = .copula_dynamic_predict(object, h = h, nsim = nsim, sim_method = sim_method[1], forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...),
                  "adcc" = .copula_dynamic_predict(object, h = h, nsim = nsim, sim_method = sim_method[1], forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...))
    return(out)
}

#' @method predict dcc.estimate
#' @rdname predict.tsmarch
#' @export
#'
#'
predict.dcc.estimate <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"),
                                    forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
{
    out <- switch(object$spec$dynamics$model,
                  "constant" = .dcc_constant_predict(object, h = h, nsim = nsim, sim_method = sim_method[1], forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...),
                  "dcc" = .dcc_dynamic_predict(object, h = h, nsim = nsim, sim_method = sim_method[1], forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...),
                  "adcc" = .dcc_dynamic_predict(object, h = h, nsim = nsim, sim_method = sim_method[1], forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...))
    return(out)
}

#' @method predict gogarch.estimate
#' @rdname predict.tsmarch
#' @export
#'
predict.gogarch.estimate <- function(object, h = 1, nsim = 1000, sim_method = c("parametric","bootstrap"),
                                     forc_dates = NULL, cond_mean = NULL, seed = NULL, ...)
{
    .gogarch_predict(object = object, h = h ,nsim = nsim, sim_method = sim_method[1],
                     forc_dates = forc_dates, cond_mean = cond_mean, seed = seed, ...)
}

# sigma ---------------------------------------------------
sigma_extractor <- function(object)
{
    if (any(grepl("estimate", class(object)))) {
        return(sigma(object$univariate))
    } else if (any(grepl("predict", class(object)))) {

    } else if (any(grepl("simulate", class(object)))) {

    } else {
        stop("\nunknown class")
    }
}

# residuals ---------------------------------------------------

standard_residuals_estimate <- function(object, standardize = FALSE, ...)
{
    out <- residuals(object$spec$univariate, standardize = standardize)
    colnames(out) <- object$spec$series_names
    return(out)
}

model_residuals_estimate <- function(object, ...)
{
    if (is(object, "cgarch.estimate")) {
        out <- object$copula_residuals
        colnames(out) <- object$spec$series_names
        out <- xts(out, object$spec$target$index)
        return(out)
    } else {
        return(standard_residuals_estimate(object = object))
    }
}

whitened_residuals_estimate <- function(object, ...)
{
    out <- object$whitened_residuals
    colnames(out) <- object$spec$series_names
    out <- xts(out, object$spec$target$index)
    return(out)
}


#' Extract Model Residuals
#'
#' @description Extract the residuals of the estimated model.
#' @param object an object of class \dQuote{cgarch.estimate},
#' \dQuote{dcc.estimate} or \dQuote{gogarch.estimate}
#' @param standardize logical. Whether to standardize the residuals by the
#' conditional volatility (only valid for the \dQuote{standard} type).
#' @param type a choice of \dQuote{standard} (default), \dQuote{model}, or
#' \dQuote{whitened} residuals. The first choice is the default and represents
#' the residuals from the first stage univariate GARCH models. The second choice
#' is only useful for the copula model and represents the residuals from the copula
#' after the transformation. In the case of the DCC model this will return the
#' standard type residuals (since they are the same). The last choice represents the
#' whitened (ZCA based) residuals which are the standard residuals multiplied by the inverse
#' of the square root of the conditional covariance matrix.
#' @note
#' In the case of the GOGARCH model, the residuals are calculated as
#' \eqn{\varepsilon  A}, where A is the mixing matrix applied
#' to the independent component residuals. These will be equal to the residuals
#' of the original series only if there is no dimensionality reduction.
#' @param ... not currently used.
#' @returns An xts matrix.
#' @aliases residuals
#' @method residuals cgarch.estimate
#' @rdname residuals.tsmarch
#' @export
#'
#'
residuals.cgarch.estimate <- function(object, standardize = FALSE, type = "standard", ...)
{
    out <- switch(type,
                  "standard" = standard_residuals_estimate(object, standardize = standardize),
                  "model" = model_residuals_estimate(object),
                  "whitened" = whitened_residuals_estimate(object))
    return(out)
}

#' @method residuals dcc.estimate
#' @rdname residuals.tsmarch
#' @export
#'
#'
residuals.dcc.estimate <- function(object, standardize = FALSE, type = "standard", ...)
{
    out <- switch(type,
                  "standard" = standard_residuals_estimate(object, standardize = standardize),
                  "model" = standard_residuals_estimate(object),
                  "whitened" = whitened_residuals_estimate(object))
    return(out)
}

#' @method residuals gogarch.estimate
#' @rdname residuals.tsmarch
#' @export
#'
residuals.gogarch.estimate <- function(object, standardize = FALSE, type = "standard", ...)
{
    A <- object$ica$A
    res <- do.call(cbind, lapply(object$univariate, function(x) coredata(residuals(x))))
    res <- res %*% A
    type <- match.arg(type, c("standard","whitened"))
    if (type == "standard") {
        if (standardize) {
            S <- tscov(object)
            S <- do.call(rbind, lapply(1:dim(S)[3], function(i) diag(S[,,i])))
            res <- res/sqrt(S)
        }
    } else if (type == "whitened") {
        S <- tscov(object)
        r <- do.call(rbind, lapply(1:dim(S)[3], function(i){
            C <- .zca(S[,,i])
            res[i,,drop = FALSE] %*% C
        }))
        res <- r
    }
    res <- xts(res, object$spec$target$index)
    colnames(res) <- object$spec$series_names
    return(res)
}


# fitted ---------------------------------------------------

fitted_estimate <- function(object, ...)
{
    out <- xts(object$mu, object$spec$target$index)
    colnames(out) <- object$spec$series_names
    return(out)
}

fitted_simulate <- function(object, ...)
{
    out <- object$mu
    attr(out,"series") <- object$series_names
    attr(out,"horizon") <- object$h
    attr(out,"draws") <- object$nsim
    return(out)
}

fitted_predict <- function(object, ...)
{
    out <- object$mu
    attr(out,"series") <- object$series_names
    attr(out,"h") <- object$h
    attr(out,"draws") <- object$nsim
    attr(out,"index")  <- as.character(object$forc_dates)
    return(out)
}


#' Extract Model Fitted Values
#'
#' @description Extract the fitted values from an object.
#' @param object an object from one of the estimated, simulated or predicted
#' model classes in the package.
#' @param ... not currently used.
#' @returns For an estimated object an \dQuote{xts} matrix of the fitted values (either
#' zero or a constant), else for simulated or predicted objects a horizon x n_series x n_draws
#' array of the distributional values. The array output includes attributes such as the
#' date index (for the predicted object), horizon, n_draws and series names.
#' @aliases fitted
#' @method fitted cgarch.estimate
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.cgarch.estimate <- function(object, ...)
{
    return(fitted_estimate(object, ...))
}

#' @aliases fitted
#' @method fitted cgarch.simulate
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.cgarch.simulate <- function(object, ...)
{
    return(fitted_simulate(object, ...))
}

#' @aliases fitted
#' @method fitted cgarch.predict
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.cgarch.predict <- function(object, ...)
{
    return(fitted_predict(object, ...))
}

#' @method fitted dcc.estimate
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.dcc.estimate <- function(object, ...)
{
    return(fitted_estimate(object, ...))
}

#' @method fitted dcc.simulate
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.dcc.simulate <- function(object, ...)
{
    return(fitted_simulate(object, ...))
}

#' @method fitted dcc.predict
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.dcc.predict <- function(object, ...)
{
    return(fitted_predict(object, ...))
}



#' @method fitted gogarch.estimate
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.gogarch.estimate <- function(object, ...)
{
    return(.gogarch_fitted(object, ...))
}

# check to use gogarch_fitted instead

#' @method fitted gogarch.predict
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.gogarch.predict <- function(object, ...)
{
    return(fitted_predict(object, ...))
}

#' @method fitted gogarch.simulate
#' @rdname fitted.tsmarch
#' @export
#'
#'
fitted.gogarch.simulate <- function(object, ...)
{
    return(fitted_predict(object, ...))
}

# tscov ---------------------------------------------------

tscov_dcc_estimate <- function(object, ...)
{
    n_series <- .lower_tri_dimension(length(object$H[1,]), diag = TRUE)
    H <- .tril2sym(object$H, n_series, TRUE)
    attr(H, "index") <- as.character(object$spec$target$index)
    attr(H, "series") <- object$spec$series_names
    return(H)
}

tscov_dcc_simulate <- function(object, distribution = TRUE, ...)
{
    n_series <- object$n_series
    h <- object$h
    nsim <- object$nsim
    dims <- c(n_series, n_series, h, nsim)
    H <- array(0, dim = dims)
    # take this down to C++
    for (i in 1:nsim) {
        H[,,,i] <- .tril2sym(.retain_dimensions_array(object$H, i), n_series, TRUE)
    }
    if (!distribution) {
        H <- array_matrix_mean(H)
    }
    attr(H, "index") <- 1:nsim
    attr(H, "series") <- object$series_names

    return(H)
}

tscov_dcc_predict <- function(object, distribution = TRUE, ...)
{
    n_series <- object$n_series
    h <- object$h
    nsim <- object$nsim
    dims <- c(n_series, n_series, h, nsim)
    H <- array(0, dim = dims)
    # take this down to C++
    for (i in 1:nsim) {
        H[,,,i] <- .tril2sym(.retain_dimensions_array(object$H, i), n_series, TRUE)
    }
    if (!distribution) {
        H <- array_matrix_mean(H)
    }
    attr(H,"index")  <- as.character(object$forc_dates)
    attr(H,"series") <- object$series_names
    return(H)
}

tscov_gogarch_estimate <- function(object, distribution = FALSE)
{
    A <- object$ica$A
    n_series <- .lower_tri_dimension(length(object$H[1,]), diag = TRUE)
    H <- .tril2sym(object$H, n_series, TRUE)
    attr(H,"index") <- object$spec$target$index
    attr(H,"series") <- object$spec$series_names
    return(H)
}

tscov_gogarch_predict <- function(object, distribution = FALSE)
{
    nsim <- object$nsim
    n_series <- object$spec$n_series
    h <- object$h
    H <- array(0, dim = c(n_series, n_series, h, nsim))
    for (i in 1:nsim) {
        V <-  do.call(cbind, lapply(object$univariate, function(x) x$sigma_sim[i,]^2))
        tmp <- .gogarch_covariance(V, object$ica$A)
        tmp <- .tril2sym(tmp, object$spec$n_series, TRUE)
        H[,,,i] <- tmp
    }
    if (!distribution) {
        V <- do.call(cbind, lapply(object$univariate, function(x) coredata(x$sigma^2)))
        tmp <- .gogarch_covariance(V, object$ica$A)
        H <- .tril2sym(tmp, object$spec$n_series, TRUE)
    }
    attr(H,"index") <- as.character(object$forc_dates)
    attr(H,"series") <- object$spec$series_names
    return(H)
}

tscov_gogarch_simulate <- function(object, distribution = FALSE)
{
    nsim <- object$nsim
    n_series <- object$spec$n_series
    h <- object$h
    H <- array(0, dim = c(n_series, n_series, h, nsim))
    for (i in 1:nsim) {
        V <-  do.call(cbind, lapply(object$univariate, function(x) x$sigma[i,]^2))
        tmp <- .gogarch_covariance(V, object$ica$A)
        tmp <- .tril2sym(tmp, object$spec$n_series, TRUE)
        H[,,,i] <- tmp
    }
    if (!distribution) {
        H <- array_matrix_mean(H)
    }
    attr(H,"index") <- 1:nsim
    attr(H,"series") <- object$spec$series_names
    return(H)
}


#' Covariance Extractor
#' @description
#' Extracts the conditional covariance matrices.
#' @param object an object class from one of the models in the package.
#' @param distribution whether to return the full simulated covariance distribution
#' for the predicted and simulated objects, else the average covariance across each
#' horizon.
#' @param ... none
#' @returns the covariance (see details).
#' @details
#' ## Estimation Object
#' An array of covariance matrices with time as the third dimension.
#' The returned object has attributes \sQuote{index} representing the datetime
#' and \sQuote{series} representing the series names.
#' ## Simulation and Prediction Objects
#' A 4-d array of dimensions (n_series x n_series x horizon x n_draws). If
#' \code{distribution} is FALSE, then the average covariance across all draws, an
#' array of dimensions (n_series x n_series x horizon).
#' @aliases tscov
#' @method tscov cgarch.estimate
#' @rdname tscov.tsmarch
#' @author Alexios Galanos
#' @export
#'
tscov.cgarch.estimate <- function(object, distribution = TRUE, ...)
{
    return(tscov_dcc_estimate(object = object, distribution = distribution, ...))
}

#' @method tscov cgarch.simulate
#' @rdname tscov.tsmarch
#' @export
#'
tscov.cgarch.simulate <- function(object, distribution = TRUE, ...)
{
    return(tscov_dcc_simulate(object = object, distribution = distribution, ...))
}

#' @method tscov cgarch.predict
#' @rdname tscov.tsmarch
#' @export
#'
tscov.cgarch.predict <- function(object, distribution = TRUE, ...)
{
    return(tscov_dcc_predict(object = object, distribution = distribution, ...))
}

#' @method tscov dcc.estimate
#' @rdname tscov.tsmarch
#' @export
#'
tscov.dcc.estimate <- function(object, distribution = TRUE, ...)
{
    return(tscov_dcc_estimate(object = object, distribution = distribution, ...))
}

#' @method tscov dcc.simulate
#' @rdname tscov.tsmarch
#' @export
#'
tscov.dcc.simulate <- function(object, distribution = TRUE, ...)
{
    return(tscov_dcc_simulate(object = object, distribution = distribution, ...))
}

#' @method tscov dcc.predict
#' @rdname tscov.tsmarch
#' @export
#'
tscov.dcc.predict <- function(object, distribution = TRUE, ...)
{
    return(tscov_dcc_predict(object = object, distribution = distribution, ...))
}

#' @method tscov gogarch.estimate
#' @rdname tscov.tsmarch
#' @export
#'
tscov.gogarch.estimate <- function(object, distribution = TRUE, ...)
{
    return(tscov_gogarch_estimate(object = object, distribution = distribution, ...))
}


#' @method tscov gogarch.predict
#' @rdname tscov.tsmarch
#' @export
#'
tscov.gogarch.predict <- function(object, distribution = TRUE, ...)
{
    return(tscov_gogarch_predict(object = object, distribution = distribution, ...))
}

#' @method tscov gogarch.simulate
#' @rdname tscov.tsmarch
#' @export
#'
tscov.gogarch.simulate <- function(object, distribution = TRUE, ...)
{
    return(tscov_gogarch_simulate(object = object, distribution = distribution, ...))
}


# tscor ---------------------------------------------------

tscor_dcc_estimate <- function(object, distribution = TRUE, ...)
{
    if (any(grepl("constant", class(object)))) {
        n_series <- .lower_tri_dimension(length(object$R), diag = FALSE)
        dims <- c(n_series, n_series)
        R <- .compose_tri_matrix(object$R, dims = dims, diag = FALSE)
        colnames(R) <- rownames(R) <- object$spec$series_names
    } else {
        n_series <- .lower_tri_dimension(length(object$R[1,]), diag = FALSE)
        n <- length(object$R[1,])
        R <- .tril2sym(object$R, n_series, FALSE)
        attr(R,"index") <- as.character(object$spec$target$index)
        attr(R,"series") <- object$spec$series_names
    }
    return(R)
}

tscor_dcc_simulate <- function(object, distribution = TRUE, ...)
{
    if (any(grepl("constant", class(object)))) {
        n_series <- .lower_tri_dimension(length(object$R), diag = FALSE)
        dims <- c(n_series, n_series)
        R <- .compose_tri_matrix(object$R, dims = dims, diag = FALSE)
        colnames(R) <- rownames(R) <- object$series_names
    } else {
        n_series <- object$n_series
        h <- object$h
        nsim <- object$nsim
        dims <- c(n_series, n_series, h, nsim)
        R <- array(0, dim = dims)
        for (i in 1:nsim) {
            R[,,,i] <- .tril2sym(.retain_dimensions_array(object$R, i), n_series, FALSE)
        }
        if (!distribution) {
            R <- array_matrix_mean(R)
        }
        attr(R, "index") <- 1:nsim
        attr(R,"series") <- object$series_names
    }
    return(R)
}

tscor_dcc_predict <- function(object, distribution = TRUE, ...)
{
    if (any(grepl("constant", class(object)))) {
        n_series <- .lower_tri_dimension(length(object$R), diag = FALSE)
        dims <- c(n_series, n_series)
        R <- .compose_tri_matrix(object$R, dims = dims, diag = FALSE)
        colnames(R) <- rownames(R) <- object$series_names
    } else {
        n_series <- object$n_series
        h <- object$h
        nsim <- object$nsim
        dims <- c(n_series, n_series, h, nsim)
        R <- array(0, dim = dims)
        for (i in 1:nsim) {
            R[,,,i] <- .tril2sym(.retain_dimensions_array(object$R, i), n_series, FALSE)
        }
        if (!distribution) {
            R <- array_matrix_mean(R)
        }
        attr(R,"series") <- object$series_names
        attr(R,"index") <- as.character(object$forc_dates)
    }
    return(R)
}


tscor_gogarch_estimate <- function(object, distribution = TRUE, ...)
{
    n_series <- .lower_tri_dimension(length(object$R[1,]), diag = FALSE)
    R <- .tril2sym(object$R, n_series, FALSE)
    attr(R,"index") <- object$spec$target$index
    attr(R,"series") <- object$spec$series_names
    return(R)
}

tscor_gogarch_predict <- function(object, distribution = TRUE, ...)
{
    A <- object$ica$A
    n_series <- object$spec$n_series
    h <- object$h
    nsim <- object$nsim
    R <- array(0, dim = c(n_series, n_series, h, nsim))
    for (i in 1:nsim) {
        V <-  do.call(cbind, lapply(object$univariate, function(x) x$sigma_sim[i,]^2))
        tmp <- .gogarch_correlation(V, A)
        tmp <- .tril2sym(tmp, n_series, FALSE)
        R[,,,i] <- tmp
    }
    if (!distribution) {
        H <- tscov_gogarch_predict(object, distribution = FALSE)
        R <- array(0, dim = c(n_series, n_series, h))
        for (i in 1:h) {
            R[,,i] <- cov2cor(H[,,i])
        }
    }
    attr(R,"index") <- as.character(object$forc_dates)
    attr(R,"series") <- object$spec$series_names
    return(R)
}

tscor_gogarch_simulate <- function(object, distribution = TRUE, ...)
{
    A <- object$ica$A
    n_series <- object$spec$n_series
    h <- object$h
    nsim <- object$nsim
    R <- array(0, dim = c(n_series, n_series, h, nsim))
    for (i in 1:nsim) {
        V <-  do.call(cbind, lapply(object$univariate, function(x) x$sigma[i,]^2))
        tmp <- .gogarch_correlation(V, A)
        tmp <- .tril2sym(tmp, n_series, FALSE)
        R[,,,i] <- tmp
    }
    if (!distribution) {
        R <- array_matrix_mean(R)
        for (i in 1:h) {
            diag(R[,,i]) <- 1
        }
    }
    attr(R,"index") <- 1:nsim
    attr(R,"series") <- object$spec$series_names
    return(R)
}


#' Correlation Extractor
#' @description
#' Extracts the conditional correlation matrices.
#' @param object an object class from one of the models in the package.
#' @param distribution whether to return the full simulated correlation distribution
#' for the predicted and simulated objects, else the average covariance across each
#' horizon.
#' @param ... none
#' @returns the correlation (see details).
#' @details
#' ## Estimation Object
#' An array of correlation matrices with time as the third dimension.
#' The returned object has attributes \sQuote{index} representing the datetime
#' and \sQuote{series} representing the series names.
#' ## Simulation and Prediction Objects
#' A 4-d array of dimensions (n_series x n_series x horizon x n_draws). If
#' \code{distribution} is FALSE, then the average covariance across all draws, an
#' array of dimensions (n_series x n_series x horizon).
#' @aliases tscor
#' @method tscor cgarch.estimate
#' @rdname tscor.tsmarch
#' @author Alexios Galanos
#' @export
#'
tscor.cgarch.estimate <- function(object, distribution = TRUE, ...)
{
    return(tscor_dcc_estimate(object = object, distribution = distribution, ...))
}


#' @method tscor cgarch.simulate
#' @rdname tscor.tsmarch
#' @export
tscor.cgarch.simulate <- function(object, distribution = TRUE, ...)
{

    return(tscor_dcc_simulate(object = object, distribution = distribution, ...))
}


#' @method tscor cgarch.predict
#' @rdname tscor.tsmarch
#' @export
tscor.cgarch.predict <- function(object, distribution = TRUE, ...)
{
    return(tscor_dcc_predict(object = object, distribution = distribution, ...))
}

#' @method tscor dcc.estimate
#' @rdname tscor.tsmarch
#' @export
#'
tscor.dcc.estimate <- function(object, distribution = TRUE, ...)
{
    return(tscor_dcc_estimate(object = object, distribution = distribution, ...))
}

#' @method tscor dcc.simulate
#' @rdname tscor.tsmarch
#' @export
tscor.dcc.simulate <- function(object, distribution = TRUE, ...)
{

    return(tscor_dcc_simulate(object = object, distribution = distribution, ...))
}


#' @method tscor dcc.predict
#' @rdname tscor.tsmarch
#' @export
tscor.dcc.predict <- function(object, distribution = TRUE, ...)
{
    return(tscor_dcc_predict(object = object, distribution = distribution, ...))
}


#' @method tscor gogarch.estimate
#' @rdname tscor.tsmarch
#' @export
tscor.gogarch.estimate <- function(object, distribution = TRUE, ...)
{
    return(tscor_gogarch_estimate(object = object, distribution = distribution, ...))
}

#' @method tscor gogarch.predict
#' @rdname tscor.tsmarch
#' @export
tscor.gogarch.predict <- function(object, distribution = TRUE, ...)
{
    return(tscor_gogarch_predict(object = object, distribution = distribution, ...))
}

#' @method tscor gogarch.simulate
#' @rdname tscor.tsmarch
#' @export
tscor.gogarch.simulate <- function(object, distribution = TRUE, ...)
{
    return(tscor_gogarch_simulate(object = object, distribution = distribution, ...))
}


# tscoskew ---------------------------------------------------
#
#' Coskewness Extractor
#' @description
#' Extracts the conditional coskewness matrices.
#' @param object an object class from one of the models in the package.
#' @param index the time index (integer) from which to extract a subset of the
#' coskewness array rather than the whole time series.
#' @param distribution whether to return the full simulated coskewness distribution
#' for the predicted and simulated objects.
#' @param standardize whether to standardize the 3th co-moment so that it represents
#' the coskewness.
#' @param folded whether to return the result as a folded or unfolded array. The folded
#' array is n_series x n_series x n_series x horizon (x simulation if predicted or simulated
#' object). The unfolded array is a n_series x (n_series^2) x horizon array. Calculations
#' such as weighted co-moments are based on the unfolded array using the Kronecker operator.
#' @param ... none
#' @returns the coskewness (see details).
#' @details
#' The calculation of the coskewness array from the independent factors is very
#' expensive in terms of memory footprint as well as computation time.
#' While it does take advantage of multiple threads if required (see
#' \code{\link[RcppParallel]{setThreadOptions}}), in the case of many series this
#' will quickly become difficult for systems low RAM. Because of this, there is
#' the option to extract a specific point in time output using the \code{index} argument.
#' @aliases tscoskew
#' @method tscoskew gogarch.estimate
#' @rdname tscoskew.tsmarch
#' @author Alexios Galanos
#' @export
#'
tscoskew.gogarch.estimate <- function(object, index = NULL, distribution = FALSE, standardize = TRUE, folded = TRUE, ...)
{
    parameter <- NULL
    A <- object$ica$A
    sig <- coredata(sigma(object$univariate))
    series_names <- object$spec$series_names
    colnames(sig) <- NULL
    n_series <- object$spec$n_series
    N <- object$spec$nobs
    if (!is.null(index)) {
        index <- .check_index_estimate(object, index)
        N <- length(index)
    } else {
        index <- 1:NROW(sig)
    }
    if (object$spec$distribution != 'norm') {
        skew <- sapply(object$univariate, function(x) x$parmatrix[parameter == "skew"]$value)
        shape <- sapply(object$univariate, function(x) x$parmatrix[parameter == "shape"]$value)
        lambda <- sapply(object$univariate, function(x) x$parmatrix[parameter == "lambda"]$value)
        sk <- dskewness(object$spec$distribution, skew = skew, shape = shape, lambda = lambda)
        # convert to moments since the standardized moments do not retain their
        # geometric properties in transformation
        sk <- matrix(sk, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^3
        cs <- .gogarch_coskewness(A, sk[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
        if (folded) cs <- fold3d(cs, p = 2)
        return(cs)
    } else {
        warning("\nmultivariate normal does not have coskewness")
        return(NULL)
    }
}

#' @method tscoskew gogarch.predict
#' @rdname tscoskew.tsmarch
#' @export
#'
tscoskew.gogarch.predict <- function(object, index = NULL, distribution = FALSE, standardize = TRUE, folded = TRUE, ...)
{
    parameter <- NULL
    A <- object$ic$A
    series_names <- object$spec$series_names
    n_series <- object$spec$n_series
    nsim <- NROW(object$univariate[[1]]$sigma_sim)
    h <- NCOL(object$univariate[[1]]$sigma_sim)
    if (!is.null(index)) {
        index <- .check_index_predict(object, index)
        N <- length(index)
    } else {
        index <- 1:h
        N <- length(index)
    }
    if (object$spec$distribution != 'norm') {
        skew <- object$parmatrix[parameter == "skew"]$value
        shape <- object$parmatrix[parameter == "shape"]$value
        lambda <- object$parmatrix[parameter == "lambda"]$value
        sk <- dskewness(object$spec$distribution, skew = skew, shape = shape, lambda = lambda)
        # convert to moments since the standardized moments do not retain their
        # geometric properties in transformation
        if (!distribution) {
            cs <- array(0, dim = c(n_series, n_series * n_series, N))
            ica_factors <- length(object$univariate)
            sig <- do.call(cbind, lapply(1:ica_factors, function(j) coredata(object$univariate[[j]]$sigma)))
            sktmp <- matrix(sk, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^3
            cs <- .gogarch_coskewness(A, sktmp[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
            if (folded) cs <- fold3d(cs, p = 2)
        } else {
            cs <- array(0, dim = c(n_series, n_series * n_series, N, nsim))
            ica_factors <- length(object$univariate)
            for (i in 1:nsim) {
                sig <- do.call(cbind, lapply(1:ica_factors, function(j) object$univariate[[j]]$sigma_sim[i,]))
                sktmp <- matrix(sk, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^3
                cs[,,,i] <- .gogarch_coskewness(A, sktmp[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
            }
            if (folded) cs <- fold4d(cs, p = 2)
        }
        return(cs)
    } else {
        warning("\nmultivariate normal does not have coskewness")
        return(NULL)
    }
}

#' @method tscoskew gogarch.simulate
#' @rdname tscoskew.tsmarch
#' @export
#'
tscoskew.gogarch.simulate <- function(object, index = NULL, distribution = FALSE, standardize = TRUE, folded = TRUE, ...)
{
    parameter <- NULL
    A <- object$ic$A
    series_names <- object$spec$series_names
    n_series <- object$spec$n_series
    nsim <- NROW(object$univariate[[1]]$sigma)
    h <- NCOL(object$univariate[[1]]$sigma)
    if (!is.null(index)) {
        index <- .check_index_predict(object, index)
        N <- length(index)
    } else {
        index <- 1:h
        N <- length(index)
    }
    if (object$spec$distribution != 'norm') {
        skew <- object$parmatrix[parameter == "skew"]$value
        shape <- object$parmatrix[parameter == "shape"]$value
        lambda <- object$parmatrix[parameter == "lambda"]$value
        sk <- dskewness(object$spec$distribution, skew = skew, shape = shape, lambda = lambda)
        # convert to moments since the standardized moments do not retain their
        # geometric properties in transformation
        cs <- array(0, dim = c(n_series, n_series * n_series, N, nsim))
        ica_factors <- length(object$univariate)
        for (i in 1:nsim) {
            sig <- do.call(cbind, lapply(1:ica_factors, function(j) object$univariate[[j]]$sigma[i,]))
            sktmp <- matrix(sk, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^3
            cs[,,,i] <- .gogarch_coskewness(A, sktmp[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
        }
        if (!distribution) {
            cs <- array4d_matrix_mean(cs)
            if (folded) cs <- fold3d(cs, p = 2)
        } else {
            if (folded) cs <- fold4d(cs, p = 2)
        }
        return(cs)
    } else {
        warning("\nmultivariate normal does not have coskewness")
        return(NULL)
    }
}

# tscokurt ---------------------------------------------------

#' Cokurtosis Extractor
#' @description
#' Extracts the conditional cokurtosis matrices.
#' @param object an object class from one of the models in the package.
#' @param index the time index (integer) from which to extract a subset of the
#' cokurtosis array rather than the whole time series.
#' @param distribution whether to return the full simulated cokurtosis distribution
#' for the predicted and simulated objects, else the average cokurtosis across each
#' horizon.
#' @param standardize whether to standardize the 4th co-moment so that it represents
#' the cokurtosis.
#' @param folded whether to return the result as a folded or unfolded array. The folded
#' array is n_series x n_series x n_series x n_series x horizon (x simulation if
#' predicted or simulated object). The unfolded array is a n_series x (n_series^3) x
#' horizon array. Calculations such as weighted co-moments are based on the unfolded
#' array using the Kronecker operator.
#' @param ... none
#' @returns the cokurtosis (see details).
#' @details
#' The calculation of the cokurtosis array from the independent factors is very
#' expensive in terms of memory footprint as well as computation time.
#' While it does take advantage of multiple threads if required (see
#' \code{\link[RcppParallel]{setThreadOptions}}), in the case of many series this
#' will quickly become difficult for systems low RAM. Because of this, there is
#' the option to extract a specific point in time output using the \code{index}
#' argument.
#' @aliases tscokurt
#' @method tscokurt gogarch.estimate
#' @rdname tscokurt.tsmarch
#' @author Alexios Galanos
#' @export
#'
tscokurt.gogarch.estimate <- function(object, index = NULL, distribution = FALSE, standardize = TRUE, folded = TRUE, ...)
{
    parameter <- NULL
    A <- object$ic$A
    sig <- coredata(sigma(object$univariate))
    series_names <- object$spec$series_names
    colnames(sig) <- NULL
    n_series <- object$spec$n_series
    N <- object$spec$nobs
    if (!is.null(index)) {
        index <- .check_index_estimate(object, index)
        N <- length(index)
    } else {
        index <- 1:NROW(sig)
    }
    if (object$spec$distribution != 'norm') {
        skew <- object$parmatrix[parameter == "skew"]$value
        shape <- object$parmatrix[parameter == "shape"]$value
        lambda <- object$parmatrix[parameter == "lambda"]$value
        ku <- dkurtosis(object$spec$distribution, skew = skew, shape = shape, lambda = lambda) + 3
        # convert to moments since the standardized moments do not retain their
        # geometric properties in transformation
        ku <- matrix(ku, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^4
        ks <- .gogarch_cokurtosis(A, K = ku[index,,drop = FALSE], V = sig[index,,drop = FALSE]^2, standardize)
        if (folded) ks <- fold3d(ks, p = 3)
        return(ks)
    } else {
        warning("\nmultivariate normal does not have a cokurtosis")
        return(NULL)
    }
}

#' @method tscokurt gogarch.predict
#' @rdname tscokurt.tsmarch
#' @export
#'
tscokurt.gogarch.predict <- function(object, index = NULL, distribution = FALSE, standardize = TRUE, folded = TRUE, ...)
{
    if (object$spec$distribution == 'norm') {
        warning("\nmultivariate normal does not have a cokurtosis")
        return(NULL)
    }
    parameter <- NULL
    A <- object$ic$A
    series_names <- object$spec$series_names
    n_series <- object$spec$n_series
    nsim <- NROW(object$univariate[[1]]$sigma_sim)
    h <- NCOL(object$univariate[[1]]$sigma_sim)
    if (!is.null(index)) {
        index <- .check_index_predict(object, index)
        N <- length(index)
    } else {
        index <- 1:h
        N <- length(index)
    }
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    ku <- dkurtosis(object$spec$distribution, skew = skew, shape = shape, lambda = lambda) + 3
    ica_factors <- length(object$univariate)
    if (!distribution) {
        sig <- do.call(cbind, lapply(1:ica_factors, function(j) coredata(object$univariate[[j]]$sigma)))
        kutmp <- matrix(ku, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^4
        ks <- .gogarch_cokurtosis(A, kutmp[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
        if (folded) ks <- fold3d(ks, p = 3)
    } else {
        ks <- array(0, dim = c(n_series, n_series * n_series * n_series, N, nsim))
        for (i in 1:nsim) {
            sig <- do.call(cbind, lapply(1:ica_factors, function(j) object$univariate[[j]]$sigma_sim[i,]))
            kutmp <- matrix(ku, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^4
            ks[,,,i] <- .gogarch_cokurtosis(A, kutmp[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
        }
        if (folded) ks <- fold4d(ks, p = 3)
    }
    return(ks)
}

#' @method tscokurt gogarch.simulate
#' @rdname tscokurt.tsmarch
#' @export
#'
tscokurt.gogarch.simulate <- function(object, index = NULL, distribution = FALSE, standardize = TRUE, folded = TRUE, ...)
{
    parameter <- NULL
    A <- object$ic$A
    series_names <- object$spec$series_names
    n_series <- object$spec$n_series
    nsim <- NROW(object$univariate[[1]]$sigma)
    h <- NCOL(object$univariate[[1]]$sigma)
    if (!is.null(index)) {
        index <- .check_index_predict(object, index)
        N <- length(index)
    } else {
        index <- 1:h
        N <- length(index)
    }
    if (object$spec$distribution != 'norm') {
        skew <- object$parmatrix[parameter == "skew"]$value
        shape <- object$parmatrix[parameter == "shape"]$value
        lambda <- object$parmatrix[parameter == "lambda"]$value
        ku <- dkurtosis(object$spec$distribution, skew = skew, shape = shape, lambda = lambda) + 3
        ks <- array(0, dim = c(n_series, n_series * n_series * n_series, N, nsim))
        ica_factors <- length(object$univariate)
        for (i in 1:nsim) {
            sig <- do.call(cbind, lapply(1:ica_factors, function(j) object$univariate[[j]]$sigma[i,]))
            kutmp <- matrix(ku, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE) * sig^4
            ks[,,,i] <- .gogarch_cokurtosis(A, kutmp[index, , drop = FALSE], sig[index, , drop = FALSE]^2, standardize)
        }
        if (!distribution) {
            ks <- array4d_matrix_mean(ks)
            if (folded) ks <- fold3d(ks, p = 3)
        } else {
            if (folded) ks <- fold4d(ks, p = 3)
        }
        return(ks)
    } else {
        warning("\nmultivariate normal does not have a cokurtosis")
        return(NULL)
    }
}

# coef ---------------------------------------------------

coef_estimate <- function(object, ...)
{
    estimate <- NULL
    if (sum(object$parmatrix$estimate) == 0) return(NULL)
    cf <- object$parmatrix[estimate == 1]$value
    names(cf) <- object$parmatrix[estimate == 1]$parameter
    return(cf)
}

#' Extract Model Coefficients
#'
#' @description
#' Extract the estimated coefficients of a model.
#'
#' @param object an estimated object from one of the models in the package.
#' @param ... none
#' @returns A numeric named vector of estimated coefficients.
#' @aliases coef
#' @method coef cgarch.estimate
#' @rdname coef.tsmarch
#' @author Alexios Galanos
#' @export
#'
coef.cgarch.estimate <- function(object, ...)
{
    return(coef_estimate(object = object, ...))
}

#' @method coef dcc.estimate
#' @rdname coef.tsmarch
#' @export
#'
coef.dcc.estimate <- function(object, ...)
{
    return(coef_estimate(object = object, ...))
}

#' @method coef gogarch.estimate
#' @rdname coef.tsmarch
#' @export
#'
coef.gogarch.estimate <- function(object, ...)
{
    series <- parameter <- value <- estimate <- NULL
    if (sum(object$parmatrix$estimate) == 0) return(NULL)
    cf <- object$parmatrix[estimate == 1, list(series, parameter, value)]
    return(cf)
}

# logLik ---------------------------------------------------

logLik_estimate <- function(object, ...)
{
    out <- -1.0 * object$loglik
    attr(out, "nobs") <- object$spec$nobs
    attr(out, "df") <- object$spec$npars
    class(out) <- "logLik"
    return(out)

}


#' Extract Log-Likelihood
#'
#' @description
#' Extract the log likelihood of the model at the estimated optimal parameter values.
#'
#' @param object an estimated object from one of the models in the package.
#' @param ... none
#' @returns An object of class \dQuote{logLik} with attributes for “nobs” and \dQuote{df}.
#' The latter is equal to the number of estimated parameters, whilst the former is the
#' number of timesteps (i.e. the number of observations per series).
#' @details
#' For all models in the package, the log-likelihood is a combination of the univariate
#' log-likelihood and the multivariate log-likelihood. For the GOGARCH model, being an
#' independent factor model, this is the sum of the univariate GARCH log-likelihoods
#' plus a term for the mixing matrix. The number of parameters in the GOGARCH model
#' reported (\dQuote{df}) represents the univariate independent factor parameters
#' plus the number of parameters in the rotation matrix U of the ICA algorithm.
#' @aliases logLik
#' @method logLik cgarch.estimate
#' @rdname logLik.tsmarch
#' @author Alexios Galanos
#' @export
#'
logLik.cgarch.estimate <- function(object, ...)
{
    return(logLik_estimate(object = object, ...))
}

#' @method logLik dcc.estimate
#' @rdname logLik.tsmarch
#' @export
#'
logLik.dcc.estimate <- function(object, ...)
{
    return(logLik_estimate(object = object, ...))
}


#' @method logLik gogarch.estimate
#' @rdname logLik.tsmarch
#' @export
#'
logLik.gogarch.estimate <- function(object, ...)
{
    return(logLik_estimate(object = object, ...))
}


# vcov ---------------------------------------------------
vcov_estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    estimate <- NULL
    type <- match.arg(type[1],c("H","OP","QMLE","NW"))
    if (is.null(object$hessian)) {
        if (type %in% c("H","QMLE","NW")) {
            warning("\nmodel was estimating without returning the hessian. Using type OP instead.")
            type <- "OP"
        }
    }
    N <- nrow(estfun(object))
    if (type == "H") {
        V <- bread(object)
    } else if (type == "QMLE") {
        bread. <- bread(object)
        meat. <- meat_tsmarch(object, adjust = adjust)
        V <- bread. %*% meat. %*% bread.
    } else if (type == "OP") {
        V <- vcovOPG(object, adjust = adjust)
    } else if (type == "NW") {
        bread. <- bread(object)
        meat. <- meatHAC_tsmarch(object, adjust = adjust, ...)
        V <- bread. %*% meat. %*% bread.
    }
    par_names <- .estimated_parameter_names(object$joint_parmatrix)
    colnames(V) <- rownames(V) <- par_names
    ndcc <- length(object$parmatrix[estimate == 1]$value)
    nall <- length(par_names)
    V <- V[(nall - ndcc + 1):nall, (nall - ndcc + 1):nall, drop = FALSE]
    return(V)
}

#' The Covariance Matrix of the Estimated Parameters
#'
#' @param object an object of class \dQuote{cgarch.estimate} or \dQuote{dcc.estimate}.
#' @param adjust logical. Should a finite sample adjustment be made? This amounts
#' to multiplication with n/(n-k) where n is the number of observations and k
#' the number of estimated parameters.
#' @param type valid choices are \dQuote{H} for using the numerical hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @returns The variance-covariance matrix of the estimated parameters.
#' @method vcov cgarch.estimate
#' @aliases vcov
#' @rdname vcov.tsmarch
#' @export
#'
vcov.cgarch.estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    return(vcov_estimate(object = object, adjust = adjust, type = type, ...))
}

#' @method vcov dcc.estimate
#' @rdname vcov.tsmarch
#' @export
#'
vcov.dcc.estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    return(vcov_estimate(object = object, adjust = adjust, type = type, ...))
}


# summary ---------------------------------------------------

#' Model Estimation Summary
#'
#' @description Summary method for class \dQuote{cgarch.estimate} or \dQuote{dcc.estimate}.
#' @param object an object of class \dQuote{cgarch.estimate} or \dQuote{dcc.estimate}.
#' @param vcov_type the type of standard errors based on the vcov estimate (see \code{\link{vcov}}).
#' @param ... not currently used.
#' @return A list with summary information of class \dQuote{summary.cgarch.estimate} or
#' \dQuote{summary.dcc.estimate}.
#' @aliases summary
#' @method summary cgarch.estimate
#' @rdname summary.tsmarch
#' @export
#'
#'
summary.cgarch.estimate <- function(object, vcov_type = "OP", ...)
{
    if (object$spec$dynamics$model == "constant" & object$spec$copula == "mvn") {
        estimate <- NULL
        V <- NULL
        est <- NULL
        par_names <- NULL
        coefficients <- NULL
        n_obs <- object$spec$nobs
        n_parameters <- 0
        n_series <- object$spec$n_series
        llh <- -object$loglik
        distribution <- object$spec$copula
        transformation <- object$spec$transformation
        coefficients <- NULL
        elapsed <- object$elapsed
        out <- list(coefficients = coefficients, distribution = distribution,
                    loglikelihood = llh, n_obs = n_obs, n_parameters = n_parameters,
                    n_series = n_series,
                    AIC = AIC(object),
                    BIC = BIC(object),
                    elapsed = elapsed, conditions = NULL,
                    dynamics = toupper(object$spec$dynamics$model),
                    model = "CGARCH",
                    transform = transformation)
    } else {
        estimate <- NULL
        if (is.null(object$hessian)) {
            vcov_type <- "OP"
            warning("\nmodel estimated without return_hessian flag. Only 'OP' vcov type available.")
        }
        V <- vcov(object, type = vcov_type)
        est <- object$parmatrix[estimate == 1]$value
        m <- length(est)
        n <- NCOL(V)
        par_names <- object$parmatrix[estimate == 1]$parameter
        se <- unname(sqrt(diag(abs(V))))
        tval <- est/se
        coefficients <- cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = 2*(1 - pnorm(abs(tval))))
        rownames(coefficients) <- par_names
        n_obs <- object$spec$nobs
        n_parameters <- object$spec$npars
        n_series <- object$spec$n_series
        llh <- -object$loglik
        conditions <- object$solution$conditions[c("kkt1","kkt2","evratio")]
        distribution <- object$spec$copula
        transformation <- object$spec$transformation
        coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
        elapsed <- object$elapsed
        setnames(coefficients, "rn","term")
        out <- list(coefficients = coefficients, distribution = distribution,
                    loglikelihood = llh, n_obs = n_obs, n_parameters = n_parameters,
                    n_series = n_series,
                    AIC = AIC(object),
                    BIC = BIC(object),
                    elapsed = elapsed, conditions = conditions,
                    dynamics = toupper(object$spec$dynamics$model),
                    model = "CGARCH",
                    transform = transformation)
    }
    class(out) <- "summary.cgarch.estimate"
    return(out)
}

#' @method summary dcc.estimate
#' @rdname summary.tsmarch
#' @export
#'
#'
summary.dcc.estimate <- function(object, vcov_type = "OP", ...)
{
    if (object$spec$dynamics$model == "constant" & object$spec$distribution == "mvn") {
        estimate <- NULL
        V <- NULL
        est <- NULL
        par_names <- NULL
        coefficients <- NULL
        n_obs <- object$spec$nobs
        n_parameters <- 0
        n_series <- object$spec$n_series
        llh <- -object$loglik
        distribution <- object$spec$copula
        coefficients <- NULL
        elapsed <- object$elapsed
        out <- list(coefficients = coefficients, distribution = object$spec$distribution,
                    loglikelihood = llh, n_obs = n_obs, n_parameters = n_parameters,
                    n_series = n_series,
                    AIC = AIC(object),
                    BIC = BIC(object),
                    elapsed = elapsed, conditions = NULL,
                    dynamics = toupper(object$spec$dynamics$model),
                    model = "DCC")
    } else {
        estimate <- NULL
        if (is.null(object$hessian)) {
            vcov_type <- "OP"
            warning("\nmodel estimated without return_hessian flag. Only 'OP' vcov type available.")
        }
        V <- vcov(object, type = vcov_type)
        est <- object$parmatrix[estimate == 1]$value
        m <- length(est)
        n <- NCOL(V)
        par_names <- object$parmatrix[estimate == 1]$parameter
        se <- unname(sqrt(diag(abs(V))))
        tval <- est/se
        coefficients <- cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = 2*(1 - pnorm(abs(tval))))
        rownames(coefficients) <- par_names
        n_obs <- object$spec$nobs
        n_parameters <- object$spec$npars
        n_series <- object$spec$n_series
        llh <- -object$loglik
        conditions <- object$solution$conditions[c("kkt1","kkt2","evratio")]
        distribution <- object$spec$distribution
        coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
        elapsed <- object$elapsed
        setnames(coefficients, "rn","term")
        out <- list(coefficients = coefficients, distribution = distribution,
                    loglikelihood = llh, n_obs = n_obs, n_parameters = n_parameters,
                    n_series = n_series,
                    AIC = AIC(object),
                    BIC = BIC(object),
                    elapsed = elapsed, conditions = conditions,
                    dynamics = toupper(object$spec$dynamics$model),
                    model = "DCC")
    }
    class(out) <- "summary.dcc.estimate"
    return(out)
}

#' @method summary gogarch.estimate
#' @rdname summary.tsmarch
#' @export
#'
#'
summary.gogarch.estimate <- function(object, vcov_type = "OP", ...)
{
    term <- NULL
    s <- lapply(object$univariate, function(x) summary(x, vcov_type = vcov_type))
    coefficients <- lapply(1:length(s), function(i) {
        tmp <- s[[i]]$coefficients
        tmp[,term := paste0("[IC_",i,"]:",term)]
    })
    coefficients <- rbindlist(coefficients)
    distribution <- object$spec$distribution
    loglikelihood <- logLik(object)
    n_obs <- object$spec$nobs
    n_parameters <- object$spec$npars
    n_series <- object$spec$n_series
    factors <- length(s)
    AIC <- AIC(object)
    BIC <- BIC(object)
    elapsed <- object$elapsed
    dynamics <- object$spec$garch$model
    model <- "GOGARCH"
    U <- object$ic$U
    K <- object$ic$K
    out <- list(coefficients = coefficients, distribution = distribution,
                loglikelihood = loglikelihood, n_obs = n_obs, n_parameters = n_parameters,
                n_series = n_series, factors = factors, AIC = AIC, BIC = BIC,
                elapsed = elapsed, dynamics = dynamics, model = model, U = U, K = K)
    class(out) <- "summary.gogarch.estimate"
    return(out)
}


#' Model Estimation Summary Print method
#'
#' @description Print method for class \dQuote{summary.cgarch.estimate} or
#' \dQuote{summary.dcc.estimate}
#' @param x an object of class \dQuote{summary.cgarch.estimate} or
#' \dQuote{summary.dcc.estimate}
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param ... not currently used.
#' @return Invisibly returns the original summary object.
#' @aliases print.summary.tsmarch.estimate
#' @method print summary.cgarch.estimate
#' @rdname print.tsmarch.estimate
#' @export
#'
#'
print.summary.cgarch.estimate <- function(x, digits = max(3L, getOption("digits") - 3L),
                                           signif.stars = getOption("show.signif.stars"),
                                           ...)
{
    .print_screen_cgarch(x, digits = digits, signif.stars = signif.stars, ...)
}

#' @method print summary.dcc.estimate
#' @rdname print.tsmarch.estimate
#' @export
#'
#'
print.summary.dcc.estimate <- function(x, digits = max(3L, getOption("digits") - 3L),
                                          signif.stars = getOption("show.signif.stars"),
                                          ...)
{
    .print_screen_dcc(x, digits = digits, signif.stars = signif.stars, ...)
}

#' @method print summary.gogarch.estimate
#' @rdname print.tsmarch.estimate
#' @export
#'
#'
print.summary.gogarch.estimate <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       signif.stars = getOption("show.signif.stars"),
                                       ...)
{
    .print_screen_gogarch(x, digits = digits, signif.stars = signif.stars, ...)
}

# tsaggregate ---------------------------------------------------

#' Weighted Moments Aggregation
#'
#' @description Calculates and returns the weighted moments of the
#' estimated, simulated or predicted object.
#' @param object an object of one of the model classes in the package.
#' @param weights an optional weighting vector. If NULL will use an equal weight vector.
#' It is also possible to pass a time varying weighting matrix with time along the
#' row dimension and equal to the number of time points or horizon.
#' @param distribution for the predicted and simulated objects whether to return
#' the simulated distribution of the weighted moments else the average.
#' @param ... not currently used.
#' @returns A list with the weighted moments. For an estimated object class these
#' will be xts vectors whilst for the simulated and predicted class these will be
#' objects of class \dQuote{tsmodel.distribution} capturing the distributional
#' uncertainty and for which a default plot method exists, unless argument
#' \dQuote{distribution} is set to FALSE.
#' @aliases tsaggregate
#' @method tsaggregate cgarch.estimate
#' @rdname tsaggregate.tsmarch
#' @export
#'
#'
tsaggregate.cgarch.estimate <- function(object, weights = NULL, ...)
{
    n <- object$spec$nobs
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = ncol(w), nrow = n, byrow = TRUE)
    }
    mu <- object$spec$target$mu
    mu <- rowSums(w * mu)
    mu <- xts(mu, object$spec$target$index)
    colnames(mu) <- "mu"
    sig <- .port_sigma(object, weights)
    sig <- xts(as.numeric(sig), object$spec$target$index)
    colnames(sig) <- "sigma"
    out <- list(mu = mu, sigma = sig)
    return(out)
}

#' @method tsaggregate cgarch.simulate
#' @rdname tsaggregate.tsmarch
#' @export
#'
#'
tsaggregate.cgarch.simulate <- function(object, weights = NULL, distribution = TRUE, ...)
{
    n <- object$h
    w <- .check_weights(weights, m = object$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = ncol(w), nrow = n, byrow = TRUE)
    }
    mu <- object$mu
    wmu <- .aggregate_mu(object$mu, w)
    wsigma <- .aggregate_sigma(object$H, w)
    wsigma <- sqrt(wsigma)
    if (distribution) {
        class(wsigma) <- "tsmodel.distribution"
        attr(wsigma, "date_class") <- "numeric"
        class(wmu) <- "tsmodel.distribution"
        attr(mu, "date_class") <- "numeric"
    } else {
        wsigma <- matrix(sqrt(colMeans(wsigma^2)), nrow = 1)
        wmu <- matrix(colMeans(wmu), nrow = 1)
        colnames(wsigma) <- "sigma"
        colnames(mu) <- "mu"
    }
    return(list(mu = wmu, sigma = wsigma))
}


#' @method tsaggregate cgarch.predict
#' @rdname tsaggregate.tsmarch
#' @export
#'
#'
tsaggregate.cgarch.predict <- function(object, weights = NULL, distribution = TRUE, ...)
{
    n <- object$h
    w <- .check_weights(weights, m = object$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = ncol(w), nrow = n, byrow = TRUE)
    }
    mu <- object$mu
    wmu <- .aggregate_mu(mu, w)
    wsigma <- .aggregate_sigma(object$H, w)
    wsigma <- sqrt(wsigma)
    if (distribution) {
        class(wmu) <- "tsmodel.distribution"
        colnames(wmu) <- as.character(object$forc_dates)
        colnames(wsigma) <- as.character(object$forc_dates)
        class(wsigma) <- "tsmodel.distribution"
    } else {
        wsigma <- xts(sqrt(colMeans(wsigma^2)), object$forc_dates)
        wmu <- xts(colMeans(wmu), object$forc_dates)
        colnames(wsigma) <- "sigma"
        colnames(mu) <- "mu"
    }
    return(list(mu = wmu, sigma = wsigma))
}


#' @method tsaggregate dcc.estimate
#' @rdname tsaggregate.tsmarch
#' @export
#'
#'
tsaggregate.dcc.estimate <- function(object, weights = NULL, ...)
{
    return(tsaggregate.cgarch.estimate(object = object, weights = weights, ...))
}


#' @method tsaggregate dcc.simulate
#' @rdname tsaggregate.tsmarch
#' @export
#'
#'
tsaggregate.dcc.simulate <- function(object, weights = NULL, distribution = TRUE, ...)
{
    return(tsaggregate.cgarch.simulate(object, weights = weights, distribution = distribution, ...))
}


#' @method tsaggregate dcc.predict
#' @rdname tsaggregate.tsmarch
#' @export
#'
#'
tsaggregate.dcc.predict <- function(object, weights = NULL, distribution = TRUE, ...)
{
    return(tsaggregate.cgarch.predict(object, weights = weights, distribution = distribution, ...))
}

#' @method tsaggregate gogarch.estimate
#' @rdname tsaggregate.tsmarch
#' @export
#'
tsaggregate.gogarch.estimate <- function(object, weights = NULL, ...)
{
    n <- object$spec$nobs
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    mu <- object$mu
    wmu <- mu %*% weights
    wsigma <- .port_sigma(object, w)
    if (object$spec$distribution == 'norm') {
        moments <- matrix(0, ncol = 4, nrow = n)
        colnames(moments) <- c("mu","sigma","skewness","kurtosis")
        moments[,1] <- as.numeric(wmu)
        moments[,2] <- wsigma
    } else {
        moments <- matrix(0, ncol = 4, nrow = n)
        colnames(moments) <- c("mu","sigma","skewness","kurtosis")
        moments[,1] <- as.numeric(wmu)
        moments[,2] <- wsigma
        moments[,3] <- .gogarch_port_skewness_estimate(object, weights = w, sigma = wsigma)
        moments[,4] <- .gogarch_port_kurtosis_estimate(object, weights = w, sigma = wsigma)
    }
    moments <- xts(moments, object$spec$target$index)
    out <- list(mu = moments[,1], sigma = moments[,2], skewness = moments[,3], kurtosis = moments[,4])
    return(out)
}

#' @method tsaggregate gogarch.predict
#' @rdname tsaggregate.tsmarch
#' @export
#'
tsaggregate.gogarch.predict <- function(object, weights = NULL, distribution = TRUE, ...)
{
    n <- object$h
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    if (distribution) {
        nsim <- NROW(object$univariate[[1]]$sigma_sim)
        mu <- object$mu
        w_mu <- matrix(NA, ncol = n, nrow = nsim)
        for (i in 1:n) {
            w_mu[,i] <- t(mu[i,,]) %*% w[i,]
        }
        w_sigma <- .gogarch_port_sigma_simulate(object, w)
        forc_dates <- as.character(object$forc_dates)
        if (object$spec$distribution == 'norm') {
            w_mu <- .set_distribution_class(w_mu, forc_dates)
            w_sigma <- .set_distribution_class(w_sigma, forc_dates)
            L <- list(mu = w_mu, sigma = w_sigma)
        } else {
            w_skew <- .gogarch_port_skewness_simulate(object, w, w_sigma)
            w_kurt <- .gogarch_port_kurtosis_simulate(object, w, w_sigma)
            w_mu <- .set_distribution_class(w_mu, forc_dates)
            w_sigma <- .set_distribution_class(w_sigma, forc_dates)
            w_skew <- .set_distribution_class(w_skew, forc_dates)
            w_kurt <- .set_distribution_class(w_kurt, forc_dates)
            L <- list(mu = w_mu, sigma = w_sigma, skewness = w_skew, kurtosis = w_kurt)
        }
    } else {
        mu <- t(apply(object$mu, 1, rowMeans))
        w_mu <- matrix(NA, ncol = 1, nrow = n)
        for (i in 1:n) {
            w_mu[i,] <- sum(mu[i,] * w[i,])
        }
        H <- tscov(object, distribution = FALSE)
        w_sigma <- matrix(NA, ncol = 1, nrow = n)
        for (i in 1:n) {
            w_sigma[i,] <- sqrt(w[i,] %*% H[,,i] %*% w[i,])
        }
        A <- object$ica$A
        m <- NROW(A)
        sig <- do.call(cbind, lapply(1:m, function(i) coredata(object$univariate[[i]]$sigma)))
        colnames(sig) <- NULL
        sk <- .gogarch_dskewness(object)
        sk <- matrix(sk, ncol = m, nrow = n, byrow = TRUE) * sig^3
        w_skew <- .gogarch_skewness_weighted(A, sk, w)/w_sigma^3
        ku <- .gogarch_dkurtosis(object)
        ku <- matrix(ku, ncol = m, nrow = n, byrow = TRUE) * sig^4
        w_kurt <- .gogarch_kurtosis_weighted(A, K = ku, V = sig^2, w = w)/w_sigma^4
        w_mu <- .make_xts(w_mu, object$forc_dates)
        w_sigma <- .make_xts(w_sigma, object$forc_dates)
        w_skew <- .make_xts(w_skew, object$forc_dates)
        w_kurt <- .make_xts(w_kurt, object$forc_dates)
        L <- list(mu = w_mu, sigma = w_sigma, skewness = w_skew, kurtosis = w_kurt)
        # cbind |> do.call(tmp)
    }
    return(L)
}

#' @method tsaggregate gogarch.simulate
#' @rdname tsaggregate.tsmarch
#' @export
#'
tsaggregate.gogarch.simulate <- function(object, weights = NULL, distribution = TRUE, ...)
{
    n <- object$h
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    if (distribution) {
        h <- NCOL(object$univariate[[1]]$sigma)
        nsim <- NROW(object$univariate[[1]]$sigma)
        mu <- object$mu
        w_mu <- matrix(NA, ncol = n, nrow = nsim)
        for (i in 1:n) {
            w_mu[,i] <- t(mu[i,,]) %*% w[i,]
        }
        w_sigma <- .gogarch_port_sigma_simulate(object, w)
        forc_dates <- as.character(object$forc_dates)
        if (object$spec$distribution == 'norm') {
            w_mu <- .set_distribution_class(w_mu, forc_dates)
            w_sigma <- .set_distribution_class(w_sigma, forc_dates)
            L <- list(mu = w_mu, sigma = w_sigma)
        } else {
            w_skew <- .gogarch_port_skewness_simulate(object, w, w_sigma)
            w_kurt <- .gogarch_port_kurtosis_simulate(object, w, w_sigma)
            w_mu <- .set_distribution_class(w_mu, forc_dates)
            w_sigma <- .set_distribution_class(w_sigma, forc_dates)
            w_skew <- .set_distribution_class(w_skew, forc_dates)
            w_kurt <- .set_distribution_class(w_kurt, forc_dates)
            L <- list(mu = w_mu, sigma = w_sigma, skewness = w_skew, kurtosis = w_kurt)
        }
    } else {
        mu <- t(apply(object$mu, 1, rowMeans))
        w_mu <- matrix(NA, ncol = 1, nrow = n)
        for (i in 1:n) {
            w_mu[i,] <- sum(mu[i,] * w[i,])
        }
        H <- tscov(object, distribution = FALSE)
        w_sigma <- matrix(NA, ncol = 1, nrow = n)
        for (i in 1:n) {
            w_sigma[i,] <- sqrt(w[i,] %*% H[,,i] %*% w[i,])
        }
        A <- object$ica$A
        m <- NROW(A)
        sig <- do.call(cbind, lapply(1:m, function(i) sqrt(colMeans(coredata(object$univariate[[i]]$sigma^2)))))
        colnames(sig) <- NULL
        sk <- .gogarch_dskewness(object)
        sk <- matrix(sk, ncol = m, nrow = n, byrow = TRUE) * sig^3
        w_skew <- .gogarch_skewness_weighted(A, sk, w)/w_sigma^3
        ku <- .gogarch_dkurtosis(object)
        ku <- matrix(ku, ncol = m, nrow = n, byrow = TRUE) * sig^4
        w_kurt <- .gogarch_kurtosis_weighted(A, K = ku, V = sig^2, w = w)/w_sigma^4
        L <- list(mu = w_mu, sigma = w_sigma, skewness = w_skew, kurtosis = w_kurt)
    }
    return(L)
}




# newsimpact ---------------------------------------------------

#' News Impact Surface
#'
#' @description News impact surface of a model
#' @param object an object of one of the estimated model classes in the package.
#' @param epsilon not used.
#' @param pair the pair of series for which to generate the news impact surface.
#' @param factor the pair of factors for which to generate the news impact surface for the GOGARCH model.
#' @param type either \dQuote{correlation} or \dQuote{covariance}.
#' @param ... additional parameters passed to the method.
#' @returns An object of class \dQuote{tsmarch.newsimpact} which has a plot method.
#' @method newsimpact cgarch.estimate
#' @rdname newsimpact.tsmarch
#' @export
#'
#
newsimpact.cgarch.estimate <- function(object, epsilon = NULL, pair = c(1,2), type = "correlation", ...)
{
    type <- match.arg(type[1], c("correlation","covariance"))
    max_pair <- object$spec$n_series
    if (length(pair) != 2) stop("\npair must be a vector of length 2.")
    if (max(pair) > max_pair) stop(paste0("\nmax(pair) must not exceed number of series (", max_pair,")"))
    switch(type,
           "correlation" = .copula_newsimpact_correlation(object, pair),
           "covariance" = .copula_newsimpact_covariance(object, pair))
}

#' @method newsimpact dcc.estimate
#' @rdname newsimpact.tsmarch
#' @export
#'
#
newsimpact.dcc.estimate <- function(object, epsilon = NULL, pair = c(1,2), type = "correlation", ...)
{
    type <- match.arg(type[1], c("correlation","covariance"))
    max_pair <- object$spec$n_series
    if (length(pair) != 2) stop("\npair must be a vector of length 2.")
    if (max(pair) > max_pair) stop(paste0("\nmax(pair) must not exceed number of series (", max_pair,")"))
    switch(type,
           "correlation" = .dcc_newsimpact_correlation(object, pair),
           "covariance" = .dcc_newsimpact_covariance(object, pair))
}


#' @method newsimpact gogarch.estimate
#' @rdname newsimpact.tsmarch
#' @export
#'
#
newsimpact.gogarch.estimate <- function(object, epsilon = NULL, pair = c(1,2), factor = c(1,1), type = "correlation", ...)
{
    type <- match.arg(type[1], c("correlation","covariance"))
    max_pair <- object$spec$n_series
    max_factor <- NROW(object$ic$A)
    if (length(pair) != 2) stop("\npair must be a vector of length 2.")
    if (max(pair) > max_pair) stop(paste0("\nmax(pair) must not exceed number of series (", max_pair,")"))

    if (length(factor) != 2) stop("\nfactor must be a vector of length 2.")
    if (max(factor) > max_factor) stop(paste0("\nmax(factor) must not exceed number of factors (", max_factor,")"))


    switch(type,
           "correlation" = .gogarch_newsimpact_correlation(object, pair, factor),
           "covariance" = .gogarch_newsimpact_covariance(object, pair, factor))
}


# plot newsimpact ---------------------------------------------------



#' News Impact Surface Plot
#'
#' @description Plot method for newsimpact class.
#' @param x an object of class \dQuote{tsmarch.newsimpact}.
#' @param y not used.
#' @param ... not currently supported.
#' @returns a plot of the newsimpact surface
#' @method plot tsmarch.newsimpact
#' @rdname plot
#' @export
#'
#
plot.tsmarch.newsimpact <- function(x, y = NULL, ...)
{

    if (x$model == "gogarch") {
        .gogarch_news_impact_surface(x, ...)
    } else {
        .dcc_news_impact_surface(x, ...)
    }
}

# tsconvolve ---------------------------------------------------

#' Convolution
#'
#' @description Generates the weighted density of the GOGARCH NIG or GH model.
#' @param object an object of class \dQuote{gogarch.estimate}, \dQuote{gogarch.predict}
#' or \dQuote{gogarch.simulate}.
#' @param weights A vector of weights of length equal to the number of series. If
#' NULL then an equal weight vector is used. A time varying matrix of weights is
#' also allowed with the correct number of rows (time points or horizon).
#' @param fft_step determines the step size for tuning the characteristic function inversion.
#' @param fft_by determines the resolution for the equally spaced support given by \code{fft_support}.
#' @param fft_support allows either a fixed support range to be given for the inversion else this is
#' calculated (if NULL) by examining the upper and lower quantiles of each independent factor modeled.
#' For the Generalized Hyperbolic distribution, it is not recommended to leave this as NULL since
#' it is quite expensive to calculate the quantiles and will significantly slow down execution time.
#' @param distribution for the simulated and predicted object, whether to apply to each draw or on the
#' average across draws (for the predicted object this is the analytic solution rather than the
#' average).
#' @param ... not currently supported.
#' @returns an object of class \dQuote{gogarch.fft} or \dQuote{gogarch.fftsim}.
#' @details
#' The Fast Fourier Transformation (FFT) is used to approximate the weighted
#' density based on its characteristic function. The weighted density is based
#' on the convolution of the scaled densities of the independent factors, by using
#' the Jacobian transformation (for more details see the vignette).
#' The returned object will be a list with the convoluted density for each time point (or each
#' time point and draw). This can then be passed to the \code{\link{dfft}}, \code{\link{pfft}} or
#' \code{\link{qfft}} methods which create smooth distributional functions.
#' @method tsconvolve gogarch.estimate
#' @aliases tsconvolve
#' @rdname tsconvolve
#' @export
#'
tsconvolve.gogarch.estimate <- function(object, weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1), ...)
{
    n <- object$spec$nobs
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
            ghyptransform(0, 1, skew = parmatrix[1,i], shape = parmatrix[2,i], parmatrix[3,i])
        }))
        pmatrix_std <- cbind(pmatrix_std, parmatrix[3,])
        colnames(pmatrix_std)[5] <- "lambda"
        pmatrix_std <- pmatrix_std[,c(4,3,2,1,5)]

    } else {
        stop("\nconvolution not required for normal distribution.")
    }
    sig <- coredata(sigma(object$univariate))[1:n, ]
    mu <- object$mu
    A <- object$ica$A
    n <- NROW(mu)
    m <- NROW(A)
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

#' @method tsconvolve gogarch.predict
#' @rdname tsconvolve
#' @export
#'
tsconvolve.gogarch.predict <- function(object, weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1), distribution = FALSE, ...)
{
    if (distribution) {
        L <- .gogarch_convolution_prediction_d(object,  weights = weights, fft_step = fft_step, fft_by = fft_by, fft_support = fft_support)
    } else {
        L <- .gogarch_convolution_prediction(object,  weights = weights, fft_step = fft_step, fft_by = fft_by, fft_support = fft_support)
    }
    return(L)
}

#' @method tsconvolve gogarch.simulate
#' @rdname tsconvolve
#' @export
#'
tsconvolve.gogarch.simulate <- function(object, weights = NULL, fft_step = 0.001, fft_by = 0.0001, fft_support = c(-1, 1), distribution = FALSE, ...)
{
    if (distribution) {
        L <- .gogarch_convolution_simulation_d(object,  weights = weights, fft_step = fft_step, fft_by = fft_by, fft_support = fft_support)
    } else {
        L <- .gogarch_convolution_simulation(object,  weights = weights, fft_step = fft_step, fft_by = fft_by, fft_support = fft_support)
    }
    return(L)
}

# fft density ---------------------------------------------------

#' FFT density, distribution and quantile method
#'
#' @param object an object of class \dQuote{gogarch.fft} formed by calling
#' the \code{\link{tsconvolve}} method.
#' @param index the time index on which to generate the function.
#' @param sim the simulation draw on which to generate the function (for the
#' predict and simulate objects).
#' @param ... additional parameters passed to the method.
#' @details These methods generate smooth approximation to the distribution
#' functions, returning in turn a function object which can be further called
#' with additional arguments expected from the density, distribution and quantile
#' functions. Random number generation can be achieved through the use of uniform
#' random variates and the quantile function (inverse CDF method).
#' @return an function object which can then be called to calculate the density
#' distribution or quantile for a particular instance (or time point),
#' @export
#' @rdname dist_fft
#'
dfft <- function(object, ...)
{
    UseMethod("dfft")
}

#' @export
#' @rdname dist_fft
#'
pfft <- function(object, ...)
{
    UseMethod("pfft")
}

#' @export
#' @rdname dist_fft
#'
qfft <- function(object, ...)
{
    UseMethod("qfft")
}

#' @export
#' @method dfft gogarch.fft
#' @rdname dist_fft
#'
dfft.gogarch.fft <- function(object, index = 1, ...)
{
    y <- object$y[[index]]
    support <- attr(y, "support")
    x <- seq(support[1], support[2], by = object$fft_by)
    f <- .create_pdf(x, dx = y, h = object$fft_by, mu = object$mu[index])
    return(f)
}

#' @export
#' @method pfft gogarch.fft
#' @rdname dist_fft
#'
pfft.gogarch.fft <- function(object, index = 1, ...)
{
    y <- object$y[[index]]
    support <- attr(y, "support")
    x <- seq(support[1], support[2], by = object$fft_by)
    f <- .create_cdf(x, dx = y, h = object$fft_by, mu = object$mu[index])
    return(f)
}

#' @export
#' @method qfft gogarch.fft
#' @rdname dist_fft
#'
qfft.gogarch.fft <- function(object, index = 1, ...)
{
    y <- object$y[[index]]
    support <- attr(y, "support")
    x <- seq(support[1], support[2], by = object$fft_by)
    f <- .create_icdf(x, dx = y, h = object$fft_by, mu = object$mu[index])
    return(f)
}

#' @export
#' @method dfft gogarch.fftsim
#' @rdname dist_fft
#'
dfft.gogarch.fftsim <- function(object, index = 1, sim = 1, ...)
{
    n_draws <- length(object$y)
    n_index <- length(object$y[[1]])
    index <- as.integer(index)
    sim <- as.integer(sim)
    if (index > n_index | index <= 0) stop("\nindex out of bounds")
    if (sim > n_draws | sim <= 0) stop("\nsim out of bounds")
    y <- object$y[[sim]][[index]]
    rnge <- attr(y, "support")
    support <- seq(rnge[1], rnge[2], by = object$fft_by)
    f <- .create_pdf(x = support, dx = y, h = object$fft_by, mu = object$mu[sim, index])
    return(f)
}

#' @export
#' @method pfft gogarch.fftsim
#' @rdname dist_fft
#'
pfft.gogarch.fftsim <- function(object, index = 1, sim = 1, ...)
{
    n_draws <- length(object$y)
    n_index <- length(object$y[[1]])
    index <- as.integer(index)
    sim <- as.integer(sim)
    if (index > n_index | index <= 0) stop("\nindex out of bounds")
    if (sim > n_draws | sim <= 0) stop("\nsim out of bounds")
    y <- object$y[[sim]][[index]]
    rnge <- attr(y, "support")
    support <- seq(rnge[1], rnge[2], by = object$fft_by)
    f <- .create_cdf(x = support, dx = y, h = object$fft_by, mu = object$mu[sim, index])
    return(f)
}

#' @export
#' @method qfft gogarch.fftsim
#' @rdname dist_fft
#'
qfft.gogarch.fftsim <- function(object, index = 1, sim = 1, ...)
{
    n_draws <- length(object$y)
    n_index <- length(object$y[[1]])
    index <- as.integer(index)
    sim <- as.integer(sim)
    if (index > n_index | index <= 0) stop("\nindex out of bounds")
    if (sim > n_draws | sim <= 0) stop("\nsim out of bounds")
    y <- object$y[[sim]][[index]]
    rnge <- attr(y, "support")
    support <- seq(rnge[1], rnge[2], by = object$fft_by)
    f <- .create_icdf(x = support, dx = y, h = object$fft_by, mu = object$mu[sim, index])
    return(f)
}

# plots ---------------------------------------------------

#' Dynamic Correlation Model Plots
#'
#' @param x an object of class \dQuote{dcc.estimate} or \dQuote{cgarch.estimate}.
#' @param y not used
#' @param series for the constant correlation a vector of length 2 indicating the
#' series numbers to use for the pairwise correlation plot. For the dynamic
#' correlation model, if NULL will include all series, else a vector of integers
#' of the series to include.
#' @param index_format for the dynamic correlation plot the x-axis label formatting.
#' @param labels whether to include y-axis labels. For a large number of series
#' it is best to leave this as FALSE.
#' @param cex_labels the shrink factor for the y-axis labels if included.
#' @param ... additional parameters passed to the plotting functions.
#'
#' @return plots of the correlation.
#' @export
#' @method plot dcc.estimate
#' @rdname correlation_plots
plot.dcc.estimate <- function(x, y = NULL, series = NULL, index_format = "%Y",
                              labels = FALSE, cex_labels = 0.8, ...)
{
    is_constant <- ifelse(x$spec$dynamics$model == "constant", TRUE, FALSE)
    if (is_constant) {
        .constant_bivariate_correlation_plot(x, series = series, ...)
    } else {
        .dynamic_pairwise_correlation_plot(x, series = series,
                                           index_format = index_format,
                                           labels = labels,
                                           cex_labels = cex_labels, ...)

    }
}

#' @export
#' @method plot cgarch.estimate
#' @rdname correlation_plots
plot.cgarch.estimate <- function(x, y = NULL, series = NULL, index_format = "%Y",
                              labels = FALSE, cex_labels = 0.8, ...)
{
    is_constant <- ifelse(x$spec$dynamics$model == "constant", TRUE, FALSE)
    if (is_constant) {
        .constant_bivariate_correlation_plot(x, series = series, ...)
    } else {
        .dynamic_pairwise_correlation_plot(x, series = series,
                                           index_format = index_format,
                                           labels = labels,
                                           cex_labels = cex_labels, ...)

    }
}

# value at risk ---------------------------------------------------

value_at_risk_dcc <- function(object, weights = NULL, alpha = 0.05) {
    alpha <- alpha[1]
    if (alpha <= 0 | alpha >= 1) stop("\nalpha must be between 0 and 1.")
    if (is.null(weights)) {
        warning("\nweights not specified, using equal weights.")
        weights <- rep(1/object$spec$n_series, object$spec$n_series)
    }
    a <- tsaggregate(object, weights = weights)
    var_value <- apply(a$mu, 2, quantile, alpha)
    if (any(grepl(pattern = "predict", x = class(object)))) {
        var_value <- xts(as.numeric(var_value), object$forc_dates)
    } else {
        var_value <- matrix(as.numeric(var_value), ncol = 1)
    }
    colnames(var_value) <- "VaR"
    return(var_value)
}

value_at_risk_gogarch <- function(object, weights = NULL, alpha = 0.05) {
    alpha <- alpha[1]
    if (alpha <= 0 | alpha >= 1) stop("\nalpha must be between 0 and 1.")
    if (is.null(weights)) {
        warning("\nweights not specified, using equal weights.")
        weights <- rep(1/object$spec$n_series, object$spec$n_series)
    }
    nsim <- object$nsim
    mu <- object$mu
    n <- object$h
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    w_mu <- matrix(NA, ncol = n, nrow = nsim)
    for (i in 1:n) {
        w_mu[,i] <- t(mu[i,,]) %*% w[i,]
    }
    var_value <- apply(w_mu, 2, quantile, alpha)
    if (any(grepl(pattern = "predict", x = class(object)))) {
        var_value <- xts(as.numeric(var_value), object$forc_dates)
    } else {
        var_value <- matrix(as.numeric(var_value), ncol = 1)
    }
    colnames(var_value) <- "VaR"
    return(var_value)
}


#' @rdname value_at_risk
#' @export
#'
value_at_risk <- function(object, ...)
{
    UseMethod("value_at_risk")
}

#' Value at Risk (VaR) method for predicted and simulated objects
#'
#' @param object an object generated from the predict or simulate methods.
#' @param weights a vector of weights of length equal to the number of series. If
#' NULL then an equal weight vector is used.
#' @param alpha the quantile level for the value at risk.
#' @param ... not used.
#' @return a matrix of the value at risk. For predict type input objects this will be an xts
#' matrix with index the forecast dates.
#' @export
#' @method value_at_risk gogarch.predict
#' @rdname value_at_risk
value_at_risk.gogarch.predict <- function(object, weights = NULL, alpha = 0.05, ...)
{
    var_value <- value_at_risk_gogarch(object, weights = weights, alpha = alpha)
    return(var_value)
}

#' @export
#' @method value_at_risk dcc.predict
#' @rdname value_at_risk
value_at_risk.dcc.predict <- function(object, weights = NULL, alpha = 0.05, ...)
{
    var_value <- value_at_risk_dcc(object, weights = weights, alpha = alpha)
    return(var_value)
}

#' @export
#' @method value_at_risk cgarch.predict
#' @rdname value_at_risk
value_at_risk.cgarch.predict <- function(object, weights = NULL, alpha = 0.05, ...)
{
    var_value <- value_at_risk_dcc(object, weights = weights, alpha = alpha)
    return(var_value)
}


#' @export
#' @method value_at_risk gogarch.simulate
#' @rdname value_at_risk
value_at_risk.gogarch.simulate <- function(object, weights = NULL, alpha = 0.05, ...)
{
    var_value <- value_at_risk_gogarch(object, weights = weights, alpha = alpha)
    return(var_value)
}

#' @export
#' @method value_at_risk dcc.simulate
#' @rdname value_at_risk
value_at_risk.dcc.simulate <- function(object, weights = NULL, alpha = 0.05, ...)
{
    var_value <- value_at_risk_dcc(object, weights = weights, alpha = alpha)
    return(var_value)
}

#' @export
#' @method value_at_risk cgarch.simulate
#' @rdname value_at_risk
value_at_risk.cgarch.simulate <- function(object, weights = NULL, alpha = 0.05, ...)
{
    var_value <- value_at_risk_dcc(object, weights = weights, alpha = alpha)
    return(var_value)
}


# expected shortfall ---------------------------------------------------

expected_shortfall_dcc <- function(object, weights = NULL, alpha = 0.05) {
    alpha <- alpha[1]
    if (alpha <= 0 | alpha >= 1) stop("\nalpha must be between 0 and 1.")
    if (is.null(weights)) {
        warning("\nweights not specified, using equal weights.")
        weights <- rep(1/object$spec$n_series, object$spec$n_series)
    }
    a <- tsaggregate(object, weights = weights)
    es_value <- apply(a$mu, 2, .es_empirical_calculation, alpha)
    if (any(grepl(pattern = "predict", x = class(object)))) {
        es_value <- xts(as.numeric(es_value), object$forc_dates)
    } else {
        es_value <- matrix(as.numeric(es_value), ncol = 1)
    }
    colnames(es_value) <- "ES"
    return(es_value)
}

expected_shortfall_gogarch <- function(object, weights = NULL, alpha = 0.05) {
    alpha <- alpha[1]
    if (alpha <= 0 | alpha >= 1) stop("\nalpha must be between 0 and 1.")
    if (is.null(weights)) {
        warning("\nweights not specified, using equal weights.")
        weights <- rep(1/object$spec$n_series, object$spec$n_series)
    }
    nsim <- object$nsim
    mu <- object$mu
    n <- object$h
    w <- .check_weights(weights, m = object$spec$n_series, n)
    if (NROW(w) == 1) {
        w <- matrix(w, ncol = NCOL(w), nrow = n, byrow = TRUE)
    }
    w_mu <- matrix(NA, ncol = n, nrow = nsim)
    for (i in 1:n) {
        w_mu[,i] <- t(mu[i,,]) %*% w[i,]
    }
    es_value <- apply(w_mu, 2, .es_empirical_calculation, alpha)
    if (any(grepl(pattern = "predict", x = class(object)))) {
        es_value <- xts(as.numeric(es_value), object$forc_dates)
    } else {
        es_value <- matrix(as.numeric(es_value), ncol = 1)
    }
    colnames(es_value) <- "ES"
    return(es_value)
}


#' @rdname expected_shortfall
#' @export
#'
expected_shortfall <- function(object, ...)
{
    UseMethod("expected_shortfall")
}

#' Expected Shortfall (ES) method for predicted and simulated objects
#'
#' @param object an object generated from the predict or simulate methods.
#' @param weights a vector of weights of length equal to the number of series. If
#' NULL then an equal weight vector is used.
#' @param alpha the quantile level for the value at risk.
#' for the GOGARCH model.
#' @param ... not used.
#' @return a matrix of the expected shortfall. For predict type input objects this will be an xts
#' matrix with index the forecast dates.
#' @export
#' @method expected_shortfall gogarch.predict
#' @rdname expected_shortfall
expected_shortfall.gogarch.predict <- function(object, weights = NULL, alpha = 0.05, ...)
{
    es_value <- expected_shortfall_gogarch(object, weights = weights, alpha = alpha)
    return(es_value)
}

#' @export
#' @method expected_shortfall dcc.predict
#' @rdname expected_shortfall
expected_shortfall.dcc.predict <- function(object, weights = NULL, alpha = 0.05, ...)
{
    es_value <- expected_shortfall_dcc(object, weights = weights, alpha = alpha)
    return(es_value)
}

#' @export
#' @method expected_shortfall cgarch.predict
#' @rdname expected_shortfall
expected_shortfall.cgarch.predict <- function(object, weights = NULL, alpha = 0.05, ...)
{
    es_value <- expected_shortfall_dcc(object, weights = weights, alpha = alpha)
    return(es_value)
}


#' @export
#' @method expected_shortfall gogarch.simulate
#' @rdname expected_shortfall
expected_shortfall.gogarch.simulate <- function(object, weights = NULL, alpha = 0.05, ...)
{
    es_value <- expected_shortfall_gogarch(object, weights = weights, alpha = alpha)
    return(es_value)
}

#' @export
#' @method expected_shortfall dcc.simulate
#' @rdname expected_shortfall
expected_shortfall.dcc.simulate <- function(object, weights = NULL, alpha = 0.05, ...)
{
    es_value <- expected_shortfall_dcc(object, weights = weights, alpha = alpha)
    return(es_value)
}

#' @export
#' @method expected_shortfall cgarch.simulate
#' @rdname expected_shortfall
expected_shortfall.cgarch.simulate <- function(object, weights = NULL, alpha = 0.05, ...)
{
    es_value <- expected_shortfall_dcc(object, weights = weights, alpha = alpha)
    return(es_value)
}


# pit ---------------------------------------------------

#' Probability Integral Transform (PIT) for weighted FFT densities
#'
#' @param object an object of class \dQuote{gogarch.fft} or \dQuote{gogarch.fftsim}.
#' @param actual a vector of realized values representing the weighted values of the
#' series for the period under consideration.
#' @param ... not used.
#' @return a matrix.
#' @export
#' @method pit gogarch.fft
#' @rdname pit
pit.gogarch.fft <- function(object, actual, ...)
{
    n <- length(object$y)
    p <- matrix(0, ncol = 1, nrow = n)
    if (missing(actual)) stop("\nactual cannot be missing")
    actual <- as.numeric(actual)
    if (length(actual) != n) stop(paste0("\nactual must be a vector of length ", n))
    for (i in 1:n) {
        pd <- pfft(object, index = i)
        p[i,1] <- pd(actual[i])
    }
    return(p)
}

