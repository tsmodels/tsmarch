.garch_extract_index <- function(x)
{
    out <- x[[1]]$spec$target$index
    return(out)
}

.garch_initialize_states <- function(object, method)
{
    model <- object$spec$model$model
    arma <- object$spec$model$arma
    if (is.null(arma)) arma <- c(0,0)
    maxpq <- max(object$spec$model$order, arma)
    init_variance <- arch_initial <- init_residuals <- init_std_residuals <- init_permanent_component <- NULL
    series_init <- resid_init <- NULL
    if (maxpq > 0) {
        if (method == "start") {
            init_variance <- rep(object$var_initial, maxpq)
            arch_initial <-  object$arch_initial
            if (object$spec$model$model == "cgarch") init_permanent_component <- init_variance
        } else {
            init_variance <- tail(as.numeric(sigma(object)^2), maxpq)
            init_residuals <- tail(as.numeric(residuals(object)), maxpq)
            init_std_residuals <- tail(as.numeric(residuals(object, standardize = TRUE)), maxpq)
            if (object$spec$model$model == "cgarch") init_permanent_component <- tail(object$permanent_component, maxpq)
            if (sum(arma) > 0) {
                series_init <- tail(as.numeric(object$spec$target$y), maxpq)
                resid_init <- tail(as.numeric(residuals(object)), maxpq)
            }
            arch_initial <- NULL
        }
        if (object$spec$model$model == "cgarch") {
            init_variance <- cbind(init_permanent_component, init_variance)
        }
    }
    return(list(variance = init_variance, residuals = init_residuals, std_residuals = init_std_residuals, arch_initial = arch_initial,
                series_init = series_init, resid_init = resid_init))
}

.garch_simulate_model <- function(object, nsim, h, burn, innov, init_method)
{
    init_states <- .garch_initialize_states(object, init_method)
    gspec <- object$spec
    gspec$parmatrix <- copy(object$parmatrix)
    args <- list(object = gspec, h = h, nsim = nsim, burn = burn, var_init = init_states$variance,
                 innov_init = init_states$std_residuals, innov = innov, arch_init = init_states$arch_initial)
    if (!is.null(init_states$series_init)) args$series_init <- init_states$series_init
    if (!is.null(init_states$resid_init)) args$resid_init <- init_states$resid_init
    sim <- do.call(simulate, args)
    return(sim)
}


.garch_filter_model <- function(object, y)
{
    m <- NCOL(y)
    new_fit <- NULL
    new_fit <- future_lapply(1:m, function(i) {
        tsfilter(object$spec$univariate[[i]], y = y[,i])
    }, future.packages = "tsgarch", future.seed = TRUE)
    new_fit <- eval(new_fit)
    names(new_fit) <- names(object$spec$univariate)
    new_fit <- to_multi_estimate(new_fit)
    return(new_fit)
}
