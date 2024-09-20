.print_screen_cgarch <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), table.caption = paste0(toupper(x$model)," Model Summary\n"), ...)
{
    term <- NULL
    df <- x$df
    rdf <- x$n_parameters + 1
    if (!is.null(table.caption)) cat(table.caption)
    cat("Copula Dynamics: ", x$transform," | ",upper_case_i(x$distribution, 1, 1)," | ",x$dynamics, sep = "")
    if (!is.null(x$coefficients)) {
        coefs <- copy(x$coefficients)
        coef_names <- coefs$term
        coefs <- coefs[,term := NULL]
        coefs <- as.data.frame(coefs)
        rownames(coefs) <- coef_names
        cat("\nCoefficients:\n")
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    cat("\nN:", as.integer(x$n_obs), "| Series: ", as.integer(x$n_series))
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = 2 + digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = 2 + digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = 2 + digits)))
    cat("\n")
    invisible(x)
}

.print_screen_dcc <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), table.caption = paste0(toupper(x$model)," Model Summary\n"), ...)
{
    term <- NULL
    df <- x$df
    rdf <- x$n_parameters + 1
    if (!is.null(table.caption)) cat(table.caption)
    cat("DCC Dynamics: ", x$dynamics," | ",upper_case_i(x$distribution, 1, 1), sep = "")
    if (!is.null(x$coefficients)) {
        coefs <- copy(x$coefficients)
        coef_names <- coefs$term
        coefs <- coefs[,term := NULL]
        coefs <- as.data.frame(coefs)
        rownames(coefs) <- coef_names
        cat("\nCoefficients:\n")
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    cat("\nN:", as.integer(x$n_obs), "| Series: ", as.integer(x$n_series))
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = 2 + digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = 2 + digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = 2 + digits)))
    cat("\n")
    invisible(x)
}
