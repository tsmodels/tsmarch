# accuracy of the cross derivatives in the partitioned hessian: the second
# stage rows/columns are compared against a Richardson extrapolated full
# hessian of the joint nll for the DCC dcc mvn model. self contained: does
# not use any of the objects in helper-global.R
estimate <- NULL
y <- as.xts(globalindices[, 1:3])[1:1200, ]
y <- sweep(y, 2, colMeans(y), "-")
x <- lapply(seq_len(3), function(i){
    spec <- garch_modelspec(y[,i], model = "garch")
    estimate(spec, keep_tmb = TRUE)
})
names(x) <- colnames(y)
x <- to_multi_estimate(x)

mod <- dcc_modelspec(x, dynamics = "dcc", distribution = "mvn") |>
    estimate(control = list(trace = 0), return_hessian = TRUE)

n <- ncol(mod$hessian)
k <- nrow(mod$parmatrix[estimate == 1])
i1 <- 1:(n - k)
i2 <- (n - k + 1):n
x0 <- mod$joint_parmatrix[estimate == 1]$value
sp <- mod$spec
sp$parmatrix <- mod$parmatrix
sp$joint_parmatrix <- mod$joint_parmatrix

Href <- numDeriv::hessian(function(q) tsmarch:::.dcc_dynamic_joint_nll(q, sp), x0,
                          method.args = list(d = 1e-4, r = 4, eps = 1e-8))

# two step QMLE standard errors: V = D^-1 Psi'Psi D^-T, psi = g2 - C A^-1 g1,
# built once on the package hessian and once on the reference second stage rows
stage2_se <- function(C, D)
{
    H <- as.matrix(mod$hessian)
    S <- mod$scores
    A <- H[i1, i1]
    psi <- S[, i2] - S[, i1] %*% t(C %*% solve(A))
    V <- solve(D) %*% crossprod(psi) %*% t(solve(D))
    sqrt(diag(V))
}

test_that("cross derivatives match a Richardson reference hessian",{
    se_pkg <- stage2_se(as.matrix(mod$hessian)[i2, i1], as.matrix(mod$hessian)[i2, i2])
    se_ref <- stage2_se(Href[i2, i1], Href[i2, i2])
    expect_equal(unname(se_pkg), unname(se_ref), tolerance = 4e-3)
})
