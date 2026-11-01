# fixtures: 3 series, garch(1,1) first stage, hessian returned so that all
# vcov types are exercised
y <- as.xts(globalindices[, 1:3])[1:1000, ]
y <- sweep(y, 2, colMeans(y), "-")
x <- lapply(seq_len(3), function(i){
    spec <- garch_modelspec(y[,i], model = "garch")
    estimate(spec, keep_tmb = TRUE)
})
names(x) <- colnames(y)
x <- to_multi_estimate(x)

vcov_models <- list(
    "dcc dcc mvn" = dcc_modelspec(x, dynamics = "dcc", distribution = "mvn") |>
        estimate(control = list(trace = 0), return_hessian = TRUE),
    "dcc dcc mvt" = dcc_modelspec(x, dynamics = "dcc", distribution = "mvt") |>
        estimate(control = list(trace = 0), return_hessian = TRUE),
    "dcc constant mvt" = dcc_modelspec(x, dynamics = "constant", distribution = "mvt") |>
        estimate(control = list(trace = 0), return_hessian = TRUE),
    "cgarch dcc mvn" = cgarch_modelspec(x, dynamics = "dcc", copula = "mvn", transformation = "parametric") |>
        estimate(control = list(trace = 0), return_hessian = TRUE),
    "cgarch constant mvt" = cgarch_modelspec(x, dynamics = "constant", copula = "mvt", transformation = "parametric") |>
        estimate(control = list(trace = 0), return_hessian = TRUE)
)

# the minimised second stage objective and its per-observation pieces for
# each model class
.stage2_fns <- function(mod)
{
    if (inherits(mod, "dcc.estimate")) {
        if (inherits(mod, "dcc.dynamic")) {
            list(objective = function(p, spec) .dcc_dynamic_values(p, spec, type = "nll"),
                 pieces = function(p, spec) .dcc_dynamic_values(p, spec, type = "dccll_vec"))
        } else {
            list(objective = function(p, spec) .dcc_constant_values(p, spec, type = "nll"),
                 pieces = function(p, spec) .dcc_constant_values(p, spec, type = "ll_vec"))
        }
    } else {
        if (inherits(mod, "cgarch.dynamic")) {
            list(objective = function(p, spec) .copula_dynamic_values(p, spec, type = "nll"),
                 pieces = function(p, spec) .copula_dynamic_values(p, spec, type = "ll_vec"))
        } else {
            list(objective = function(p, spec) .copula_constant_values(p, spec, type = "nll"),
                 pieces = function(p, spec) .copula_constant_values(p, spec, type = "ll_vec"))
        }
    }
}

.check_vcov <- function(mod)
{
    estimate <- NULL
    fns <- .stage2_fns(mod)
    spec <- mod$spec
    n <- ncol(mod$hessian)
    k <- nrow(mod$parmatrix[estimate == 1])
    i2 <- (n - k + 1):n
    i1 <- 1:(n - k)
    p2 <- mod$parmatrix[estimate == 1]$value
    H <- as.matrix(mod$hessian)
    S <- mod$scores
    D <- H[i2, i2]
    # (a) the second stage block is positive definite and matches a high
    # accuracy numerical hessian of the minimised objective
    Dref <- numDeriv::hessian(function(p) fns$objective(p, spec), p2,
                              method.args = list(d = 1e-4, r = 4, eps = 1e-7))
    expect_true(all(eigen(D, symmetric = TRUE)$values > 0))
    expect_true(max(abs(D - Dref))/max(abs(Dref)) < 2e-2)
    # (b) hessian type variances are positive and match the reference
    VH <- vcov(mod, type = "H")
    expect_true(all(diag(VH) > 0))
    expect_equal(unname(diag(VH)), unname(diag(solve(Dref))), tolerance = 2e-2)
    # (c) second stage score columns have the same sign as the NLL gradient
    J <- numDeriv::jacobian(function(p) fns$pieces(p, spec), p2,
                            method.args = list(d = 1e-4, r = 2, eps = 1e-7))
    J <- tail(J, nrow(S))
    for (j in seq_len(k)) {
        expect_true(cor(J[,j], S[, i2[j]]) > 0.99)
    }
    # (d) the QMLE sandwich equals the two step form D^-1 Psi'Psi D^-T with
    # psi = g2 - C A^-1 g1, and both sandwich estimators are symmetric
    A <- H[i1, i1]
    C <- H[i2, i1]
    psi <- S[, i2] - S[, i1] %*% t(C %*% solve(A))
    V2 <- solve(D) %*% crossprod(psi) %*% t(solve(D))
    VQ <- vcov(mod, type = "QMLE")
    expect_equal(unname(VQ), unname(V2), tolerance = 1e-8)
    expect_true(isSymmetric(unname(VQ)))
    VN <- vcov(mod, type = "NW")
    expect_true(isSymmetric(unname(VN)))
    # (e) OP variances are positive
    expect_true(all(diag(vcov(mod, type = "OP")) > 0))
    return(invisible(NULL))
}

for (nm in names(vcov_models)) {
    test_that(paste0("vcov hessian sign and sandwich for ", nm),{
        .check_vcov(vcov_models[[nm]])
    })
}
