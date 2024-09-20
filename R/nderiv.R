
nderiv_function_wrapper <- function(func, x, lower, upper) {
    n <- length(x)
    idx <- x == 0.0
    orders_of_magnitude <- rep(0, n)
    orders_of_magnitude[idx] <- 1.0
    orders_of_magnitude[!idx] <- 10^ceiling(log10(abs(x[!idx])))

    scaled_x <- x / orders_of_magnitude
    scaled_lower <- lower / orders_of_magnitude
    scaled_upper <- upper / orders_of_magnitude
    scaled_deltas <- rep(0, length(scaled_x))
    for (i in seq_along(scaled_x)) {
        scaled_value <- scaled_x[i]
        scaled_min_value <- scaled_lower[i]
        scaled_max_value <- scaled_upper[i]

        if (scaled_value == scaled_min_value || scaled_value == scaled_max_value) {
            stop(paste0("Value for parameter number ", i, " is on the boundary"))
        }

        if (!is.nan(scaled_min_value)) {
            distance_to_min <- scaled_value - scaled_min_value
        } else {
            distance_to_min <- Inf
        }

        if (!is.nan(scaled_max_value)) {
            distance_to_max <- scaled_max_value - scaled_value
        } else {
            distance_to_max <- Inf
        }

        if (scaled_x[i] == 0.0) {
            scaled_deltas[i] <- min(1e-5, distance_to_max / 2.5, distance_to_min / 2.5)
        } else {
            scaled_deltas[i] <- min(0.003 * abs(scaled_x[i]), distance_to_max / 2.5, distance_to_min / 2.5)
        }
    }

    function_wrapper <- function(x, ...) {
        scaled_back_x <- x * orders_of_magnitude
        tryCatch(
            func(scaled_back_x, ...),
            error = function(e) {
                stop(paste0("Cannot compute Hessian, parameters out of bounds at ", scaled_back_x))
            }
        )
    }

    return(list(function_wrapper = function_wrapper, scaled_deltas = scaled_deltas, scaled_x = scaled_x,
                orders_of_magnitude = orders_of_magnitude, n = n))
}

nderiv_jacobian <- function(func, x, lower, upper, ...) {
    sol <- nderiv_function_wrapper(func, x, lower, upper)
    f <- sol$function_wrapper
    scaled_deltas <- sol$scaled_deltas
    scaled_x <- sol$scaled_x
    orders_of_magnitude <- sol$orders_of_magnitude
    jac <- jacobian(func = f, x = scaled_x, method = "Richardson", method.args = list(d = scaled_deltas), ...)
    # Now correct back the Jacobian for the scales
    jac <- sweep(jac, 2, 1/orders_of_magnitude, FUN = "*")
    return(jac)
}

nderiv_hessian <- function(func, x, lower, upper, ...) {
    sol <- nderiv_function_wrapper(func, x, lower, upper)
    f <- sol$function_wrapper
    scaled_deltas <- sol$scaled_deltas
    scaled_x <- sol$scaled_x
    orders_of_magnitude <- sol$orders_of_magnitude
    n <- sol$n
    # Compute the Hessian matrix at best_fit_values
    hess <- hessian(f, scaled_x, method = "Richardson", method.args = list(d = scaled_deltas), ...)
    # Now correct back the Hessian for the scales
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            hess[i, j] <- hess[i, j] / (orders_of_magnitude[i] * orders_of_magnitude[j])
        }
    }
    return(hess)
}
