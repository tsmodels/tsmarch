# regression test for issue #7: first stage models with unequal numbers of
# estimated parameters (variance targeting on a subset of series) broke the
# joint hessian/score calculation because sapply returned a list
set.seed(1)
n_uneq <- 800
y_uneq <- xts(matrix(rt(n_uneq * 3, 6) * 0.01, n_uneq, 3), order.by = as.Date("2000-01-01") + seq_len(n_uneq))
colnames(y_uneq) <- c("A","B","C")

x_uneq <- suppressWarnings(lapply(1:3, function(i){
    spec <- garch_modelspec(y_uneq[,i], model = "garch", order = c(1,1), constant = TRUE,
                            distribution = "std", variance_targeting = (i == 2))
    estimate(spec, keep_tmb = TRUE)
}))
names(x_uneq) <- colnames(y_uneq)
x_uneq <- to_multi_estimate(x_uneq)

check_unequal_model <- function(spec, class_name, expect_null = FALSE) {
    mod <- estimate(spec, control = list(trace = 0))
    expect_s3_class(mod, class_name)
    k <- sum(mod$joint_parmatrix$estimate == 1)
    if (expect_null) {
        # constant Gaussian branch has no second stage parameters
        expect_null(mod$hessian)
        expect_null(mod$scores)
    } else {
        expect_equal(dim(mod$hessian), c(k, k))
        expect_true(all(is.finite(mod$hessian)))
        if (!is.null(mod$scores)) expect_equal(ncol(mod$scores), k)
    }
    # second reported symptom: scores path must also work without the hessian
    mod_nh <- estimate(spec, control = list(trace = 0), return_hessian = FALSE)
    expect_s3_class(mod_nh, class_name)
    if (!expect_null && !is.null(mod_nh$scores)) expect_equal(ncol(mod_nh$scores), k)
    return(invisible(mod))
}

test_that("first stage models have unequal parameter counts",{
    counts <- sapply(x_uneq, function(f) sum(f$parmatrix$estimate == 1))
    expect_false(length(unique(counts)) == 1)
})

test_that("cgarch dcc unequal marginals",{
    check_unequal_model(cgarch_modelspec(x_uneq, dynamics = "dcc", copula = "mvt", transformation = "parametric"), "cgarch.estimate")
    check_unequal_model(cgarch_modelspec(x_uneq, dynamics = "dcc", copula = "mvn", transformation = "parametric"), "cgarch.estimate")
})

test_that("cgarch constant unequal marginals",{
    check_unequal_model(cgarch_modelspec(x_uneq, dynamics = "constant", copula = "mvt", transformation = "parametric"), "cgarch.estimate")
    check_unequal_model(cgarch_modelspec(x_uneq, dynamics = "constant", copula = "mvn", transformation = "parametric"), "cgarch.estimate", expect_null = TRUE)
})

test_that("dcc dynamic unequal marginals",{
    check_unequal_model(dcc_modelspec(x_uneq, dynamics = "dcc", distribution = "mvt"), "dcc.estimate")
    check_unequal_model(dcc_modelspec(x_uneq, dynamics = "adcc", distribution = "mvn"), "dcc.estimate")
})

test_that("dcc constant unequal marginals",{
    check_unequal_model(dcc_modelspec(x_uneq, dynamics = "constant", distribution = "mvt"), "dcc.estimate")
    check_unequal_model(dcc_modelspec(x_uneq, dynamics = "constant", distribution = "mvn"), "dcc.estimate", expect_null = TRUE)
})
