lmatrix <- function(x, lags, constant = FALSE)
{
    lags <- lags + 1
    n <- NROW(x)
    new_x <- rbind(x, matrix(0, lags, 1))
    lag_matrix <- kronecker(matrix(1, lags, 1), new_x)
    lag_matrix <- matrix(lag_matrix[1:(NROW(lag_matrix) - lags)], nrow = (n + lags - 1), ncol = lags)
    lag_matrix <- lag_matrix[lags:n,]
    y <- lag_matrix[,1]
    x <- lag_matrix[,2:lags]
    if (constant) x = cbind(matrix(1, NROW(x), 1), x)
    return(list(y = y, x = x))
}

signif_codes <- function()
{
    return(c("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))
}


pvalue_format <- function(x) {
    z <- cut(x, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))
    as.character(z)
}

#' Engle-Sheppard Constant Correlation Test
#'
#' @description The constant correlation test of Engle and Sheppard (2001).
#' @param x an xts matrix of stationary data.
#' @param order the GARCH model order.
#' @param lags the number of lags to test for.
#' @param ... none.
#' @returns An object of class \dQuote{tstest.escc} which has a print and
#' as_flextable method.
#' @details
#' The whitened residuals from a constant correlation GARCH model are used to
#' test the hypothesis that they are IID using an artificial regression of the outer
#' products of tese residuals on a constant and their lagged outer products.
#' The null hypothesis of constant correlation can then be tested using a Wald-type
#' test on the constant and coefficients of the lagged regressors (which should
#' be zero under the null), with the test statistic asymptotically distributed
#' as Chi-squared with (lags + 1) degrees of freedom. In the returned table, the
#' parameter \dQuote{J} represents the test statistic.
#' @note
#' The test makes use of parallel functionality for estimating the univariate GARCH
#' models and the user can control the number of threads used using
#' \code{\link[future]{plan}}.
#' @aliases escc_test
#' @references
#' \insertRef{Engle2001}{tsmarch}
#' @rdname escc_test
#' @export
#'
#'
escc_test <- function(x, order = c(1,1), lags = 1, constant = FALSE, ...)
{
    n <- NROW(x)
    m <- NCOL(x)
    if (is.null(colnames(x))) colnames(x) <- paste0("S_",1:n)
    garch_model <- future_lapply(1:m, function(i) {
        g_mod <- garch_modelspec(x[,i], model = "garch", order = order, constant = constant) |> estimate(keep_tmb = TRUE)
        return(g_mod)
    }, future.packages = "tsgarch")
    garch_model <- to_multi_estimate(garch_model)
    names(garch_model) <- colnames(y)
    dcc_mod <- dcc_modelspec(garch_model, dynamics = "constant", distribution = "mvn") |> estimate()
    z <- coredata(residuals(dcc_mod, standardize = TRUE, type = "whitened"))
    outer_products <- NULL
    for (i in 1:m) {
        if (i < m) {
            for (j in (i + 1):m) {
                outer_products <- cbind(outer_products, z[,i] * z[,j])
            }
        }
    }
    j <- NCOL(outer_products)
    regressors <- regressand <- NULL
    for (i in 1:j) {
        tmp <- lmatrix(outer_products[,i,drop = FALSE], lags, TRUE)
        regressors <- rbind(regressors, tmp$x)
        regressand <- c(regressand, tmp$y)
    }
    regressors <- as.matrix(regressors)
    regressand <- as.matrix(regressand)
    data_in <- cbind(regressand, regressors[,-1])
    c_names <- paste0("x(L",1:lags,")")
    colnames(data_in) <- c("y",c_names)
    mod <- lm(y~., data = as.data.frame(data_in))
    beta <- matrix(coef(mod), ncol = 1)
    X <- t(regressors) %*% regressors
    noise <- regressand - regressors %*% beta
    sig <- (t(noise) %*% noise)/(NROW(regressors) - lags - 1)
    stat <- (t(beta) %*% X %*% beta)/sqrt(sig)
    pval <- 1 - pchisq(stat, lags + 1)
    escc_table <- as.data.table(summary(mod)$coefficients, keep.rownames = TRUE)
    setnames(escc_table,"rn","parameter")
    escc_table[,parameter := c("constant", c_names)]
    # remove constant from tables
    joint_effect_table <- data.table("parameter" = "J", "Estimate" = NA, "Std. Error" = NA, "t value" = as.numeric(stat), "Pr(>|t|)" = as.numeric(pval))
    escc_table <- rbind(escc_table, joint_effect_table)
    escc_table[,signif := pvalue_format(`Pr(>|t|)`)]
    escc_table[,'Decision(5%)' := " "]
    if (escc_table[parameter == "J"]$`Pr(>|t|)` <= 0.05) {
        escc_table[parameter == "J",'Decision(5%)' := "Reject H0"]
    } else {
        escc_table[parameter == "J",'Decision(5%)' := "Fail to Reject H0"]

    }
    symbols <- c("c", paste0("x_{t-",1:lags,"}"),"J")
    out <- list(table = escc_table,
                joint_hypothesis = c("c,x(L) = 0"),
                hypothesis = "Constant Correlation",  test_type = "Wald",
                p_value = pval, distribution = "Chi-squared", symbols = symbols, test_class = "escc",
                test_name = "Engle-Sheppard Constant Correlation",
                reference = "Engle R.F., and Sheppard K. (2001) `Theoretical and Empirical Properties of Dynamic Conditional Correlation Multivariate GARCH.` Vol. 15. UCSD Working Paper.")
    class(out) <- c("tstest.escc","tstest")
    return(out)
}


#' Print method for escc test
#'
#' @description Print method for objects inheriting class \dQuote{tstest.escc}
#' @param x an object inheriting class \dQuote{tstest.escc}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, \sQuote{significance stars} are printed.
#' @param include.decision prints out whether to reject the NULL at the 5% level
#' of significance.
#' @param ... not currently used.
#' @return Invisibly returns the original object.
#' @method print tstest.escc
#' @rdname print
#' @export
#'
#'
print.tstest.escc <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                include.decision = FALSE, ...)
{
    signif <- parameter <- `Decision(5%)` <- NULL
    cat(x$test_name)
    cat("\nHypothesis(H0) : ", x$hypothesis,"\n\n")
    tab <- copy(x$table)
    if (!include.decision) {
        tab[,`Decision(5%)` := NULL]
    }
    if (!signif.stars) {
        tab[,signif := NULL]
    } else {
        setnames(tab, "signif"," ")
    }
    cnames <- tab$parameter
    tab[,parameter := NULL]
    tab <- as.data.frame(tab)
    rownames(tab) <- cnames
    print(tab, digits = digits)
    cat("\n---")
    if (signif.stars) {
        cat("\n")
        cat(signif_codes())
    }
    return(invisible(x))
}

#' @title Transform a test object into flextable
#' @description
#' Transforms a \dQuote{tstest.escc} object into a flextable
#' with options on symbolic representation.
#' @param x an object of which inherits a \dQuote{tstest.escc} class.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, \sQuote{significance stars} are printed.
#' @param include.decision prints out whether to reject the NULL at the 5% level
#' of significance.
#' @param use.symbols for tests which either have parameters for which the
#' latex symbols were included in the calling function or for which the tests
#' generate values which can be represented as latex symbols, then these will
#' be generated.
#' @param table.caption an optional string for the table caption.
#' @param footnote.reference whether to include the reference paper of the test
#' in the footnotes.
#' @param ... not currently used. The returned object can be manipulated further
#' using flextable.
#' @return A flextable object.
#' @aliases as_flextable.tstest.escc
#' @method as_flextable tstest.escc
#' @rdname as_flextable
#' @export
#'
as_flextable.tstest.escc <- function(x, digits = max(3L, getOption("digits") - 3L),
                                         signif.stars = getOption("show.signif.stars"),
                                         include.decision = FALSE,
                                         use.symbols = TRUE, table.caption = x$test_name,
                                         footnote.reference = FALSE, ...)
{
    `Decision(5%)` <- NULL
    if (is.null(table.caption)) table.caption <- x$test_name
    tab <- copy(x$table)
    cnames <- colnames(tab)
    if (!include.decision) {
        tab[,`Decision(5%)` := NULL]
        cnames <- colnames(tab)

    }
    if (!signif.stars) {
        tab[,signif := NULL]
        cnames <- colnames(tab)
    }
    tab <- as.data.frame(tab)
    out <- flextable(tab) |> set_caption(caption = table.caption) |> align(j = "parameter", align = "left")
    if (signif.stars) {
        out <- out |> align(j = "signif", align = "left") |>
            bold(j = "signif", bold = TRUE) |>
            padding(padding.left = 0, j = "signif", part  = "all") |>
            set_header_labels(signif = "", parameter = "")
        out <- out |> add_footer_lines(values = signif_codes())
    }
    out <- out |> add_footer_lines(top = FALSE, values = c(paste0("Hypothesis(H0) : ",x$hypothesis)))
    if (use.symbols) {
        sym <- x$symbols
        for (i in 1:nrow(tab)) {
            out <- compose(out, i = i, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = i,j = 1, as_equation(sym[i]))
        }
    }
    if (footnote.reference) {
        out <- out |> add_footer_lines(top = FALSE, values = paste0("Reference: ", x$reference))
    }
    out <- colformat_double(out, digits = digits) |> autofit()
    if (signif.stars) out <- out |> hline(i = 1, part = "footer")
    return(out)
}
