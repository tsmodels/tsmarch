.print_screen_cgarch <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), table.caption = paste0(toupper(x$model)," Model Summary\n"), ...)
{
    term <- NULL
    df <- x$df
    rdf <- x$n_parameters + 1
    if (!is.null(table.caption)) cat(table.caption)
    cat("Copula Dynamics: ", x$transform," | ", .dist_names(x$distribution)," | ",x$dynamics, sep = "")
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
    cat("DCC Dynamics: ", x$dynamics," | ", .dist_names(x$distribution), sep = "")
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


.print_screen_gogarch <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), table.caption = paste0(toupper(x$model)," Model Summary\n"), ...)
{
    term <- NULL
    rdf <- x$n_parameters
    if (!is.null(table.caption)) cat(table.caption)
    cat("Factor Dynamics: ", toupper(x$dynamics)," | ", .dist_names(x$distribution), sep = "")
    if (!is.null(x$coefficients)) {
        coefs <- copy(x$coefficients)
        coef_names <- coefs$term
        coefs <- coefs[,term := NULL]
        coefs <- as.data.frame(coefs)
        rownames(coefs) <- coef_names
        cat("\nCoefficients:\n")
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    if (NCOL(x$U) < 7) {
        cat("\nRotation Matrix (U):\n")
        write.table(format(x$U, justify = "right", digits = digits), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    cat("\nN:", as.integer(x$n_obs), "| Series: ", as.integer(x$n_series),"| Factors: ", as.integer(NCOL(x$U)))
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = 2 + digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = 2 + digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = 2 + digits)))
    cat("\n")
    invisible(x)
}

.dist_names <- function(x) {
    switch(x,
           "mvn" = "MVN",
           "norm" = "MVN",
           "gaussian" = "MVN",
           "normal" = "MVN",
           "mvt" = "MVT",
           "std" = "MVT",
           "student" = "MVT",
           "nig" = "MANIG",
           "manig" = "MANIG",
           "gh" = "MAGH",
           "magh" = "MAGH")
}
