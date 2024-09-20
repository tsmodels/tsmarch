.extract_pit_parametric <- function(x)
{
    out <- coredata(do.call(cbind, lapply(x, pit)))
    colnames(out) <- names(x)
    return(out)
}

.filter_pit_parametric <- function(x)
{
    .extract_pit_parametric(x)
}

.extract_pit_empirical <- function(x)
{
    z <- coredata(do.call(cbind, lapply(x, residuals, standardize = TRUE)))
    m <- dim(z)[2]
    n <- dim(z)[1]
    ures <- matrix(NA, ncol = m, nrow = n)
    # avoid future_lapply as we are calling in parallel a function which uses
    # this in the standard error calculation
    transform_model <- lapply(1:m, function(i) {
        mod <- ecdf(sort(z[,i]))
        return(mod)
    })
    for (i in 1:m) {
        ures[,i] <- transform_model[[i]](z[,i])
    }
    return(list(ures = ures, transform_model = transform_model))
}

# n_update considers the information just prior to T+1 so as to avoid
# look-ahead bias, whereas without the n_update, the most recent estimated
# object is used which is definitely older than T+1. The reason for allowing
# n_update to be used is so that the information set can occasionally be updated
# to the most recent data, but still less than T+1
.filter_pit_empirical <- function(x, transform_model, n_update = NULL)
{
    z <- coredata(do.call(cbind, lapply(x, residuals, standardize = TRUE)))
    m <- dim(z)[2]
    n <- dim(z)[1]
    ures <- matrix(NA, ncol = m, nrow = n)
    if (!is.null(n_update)) {
        transform_model <- NULL
        transform_model <- future_lapply(1:m, function(i) {
            mod <- ecdf(sort(z[1:n_update,i]))
            return(mod)
        }, future.seed = TRUE)
        transform_model <- eval(transform_model)
    }
    for (i in 1:m) {
        ures[,i] <- transform_model[[i]](z[,i])
    }
    return(list(ures = ures, transform_model = transform_model))
}

.extract_pit_spd <- function(x, ...)
{
    z <- coredata(do.call(cbind, lapply(x, residuals, standardize = TRUE)))
    m <- dim(z)[2]
    n <- dim(z)[1]
    transform_model <- lapply(1:m, function(i) {
        spec <- spd_modelspec(z[,i], ...)
        mod <- estimate(spec, method = "pwm")
        return(mod)
    })
    ures <- matrix(NA, ncol = m, nrow = n)
    for (i in 1:m) {
        ures[,i] <- pspd(q = z[,i], object = transform_model[[i]])
    }
    return(list(ures = ures, transform_model = transform_model))
}

.filter_pit_spd <- function(x, transform_model, n_update = NULL, ...)
{
    z <- coredata(do.call(cbind, lapply(x, residuals, standardize = TRUE)))
    m <- dim(z)[2]
    n <- dim(z)[1]
    ures <- matrix(NA, ncol = m, nrow = n)
    if (!is.null(n_update)) {
        transform_model <- NULL
        transform_model <- future_lapply(1:m, function(i) {
            spec <- spd_modelspec(z[1:n_update,i], ...)
            mod <- estimate(spec, method = "pwm")
            return(mod)
        }, future.packages = "tsdistributions", future.seed = TRUE)
    }
    transform_model <- eval(transform_model)
    for (i in 1:m) {
        ures[,i] <- pspd(q = z[,i], object = transform_model[[i]])
    }
    return(list(ures = ures, transform_model = transform_model))
}

copula_transformation_estimate <- function(object, transformation) {
    if (transformation == "parametric") {
        u <- .extract_pit_parametric(object)
        transform_model <- NULL
    } else {
        tmp <- switch(transformation,
                      "spd" = .extract_pit_spd(object),
                      "empirical" = .extract_pit_empirical(object))
        u <- tmp$ures
        transform_model <- tmp$transform_model
    }
    u[u < 3.330669e-16] <- 3.330669e-16
    u[u > 0.999999] <- 0.999999
    return(list(u = u, transform_model = transform_model))
}

copula_transformation_filter <- function(object, transformation, transform_model, n_update = NULL) {
    if (transformation == "parametric") {
        u <- .extract_pit_parametric(object)
        transform_model <- NULL
    } else {
        tmp <- switch(transformation,
                      "spd" = .filter_pit_spd(x = object, transform_model = transform_model, n_update = n_update),
                      "empirical" = .filter_pit_empirical(x = object, transform_model = transform_model, n_update = n_update))
        u <- tmp$ures
        transform_model <- tmp$transform_model
    }
    u[u < 3.330669e-16] <- 3.330669e-16
    u[u > 0.999999] <- 0.999999
    return(list(u = u, transform_model = transform_model))
}

