# fixtures: small and fast, built locally (helper-global.R keeps the non-vreg globals)
# series with genuine regressor effects in the variance so the effect is material
set.seed(1042)
n_obs <- 500
vr1 <- 0.5 + runif(n_obs)
vr2 <- 0.5 + runif(n_obs)
sim_fun <- function(vreg, xi) {
    e <- rnorm(n_obs + 50)
    s <- rep(1, n_obs + 50)
    v <- c(rep(0, 50), vreg)
    for (t in 2:(n_obs + 50)) s[t] <- sqrt(0.05 + 0.1 * e[t - 1]^2 + 0.75 * s[t - 1]^2 + xi * v[t])
    tail(e * s, n_obs)
}
y_vreg <- xts(cbind(sim_fun(vr1, 0.4), sim_fun(vr2, 0.3), sim_fun(rep(0, n_obs), 0)),
              order.by = y[1:n_obs, test_series] |> index())
colnames(y_vreg) <- colnames(y)[test_series]
vreg_insample <- list(xts(vr1, order.by = index(y_vreg)),
                      xts(vr2, order.by = index(y_vreg)),
                      NULL)

x_vreg <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y_vreg[,i], model = "garch", vreg = vreg_insample[[i]])
    estimate(spec, keep_tmb = TRUE)
})
names(x_vreg) <- colnames(y_vreg)
x_vreg <- to_multi_estimate(x_vreg)

dcc_vreg_spec <- dcc_modelspec(x_vreg, dynamics = "dcc", distribution = "mvn")
dcc_vreg_estimate <- estimate(dcc_vreg_spec, control = list(trace = 0), return_hessian = FALSE)
cgarch_vreg_spec <- cgarch_modelspec(x_vreg, dynamics = "constant", transformation = "parametric", copula = "mvn")
cgarch_vreg_estimate <- estimate(cgarch_vreg_spec, control = list(trace = 0))

sn <- colnames(y_vreg)
new_dates <- seq(tail(index(y), 1), length.out = 10, by = 7)
ynew <- xts(matrix(rnorm(30), ncol = 3), order.by = new_dates)
colnames(ynew) <- sn
newvreg_filter <- list(matrix(runif(10), ncol = 1), matrix(runif(10), ncol = 1), NULL)
names(newvreg_filter) <- sn
newvreg_hor <- list(matrix(runif(3), ncol = 1), matrix(runif(3), ncol = 1), NULL)
names(newvreg_hor) <- sn
newvreg_sim <- list(matrix(runif(5), ncol = 1), matrix(runif(5), ncol = 1), NULL)
names(newvreg_sim) <- sn

# second fixture: series 1 has a 2-column variance regressor, with the two
# columns on visibly different scales so that a misordered or transposed
# %*% xi multiply would give a materially different variance intercept
vr3 <- 5 * runif(n_obs)
sim_fun_mc <- function(vreg, xi) {
    e <- rnorm(n_obs + 50)
    s <- rep(1, n_obs + 50)
    v <- rbind(matrix(0, nrow = 50, ncol = NCOL(vreg)), vreg)
    for (t in 2:(n_obs + 50)) s[t] <- sqrt(0.05 + 0.1 * e[t - 1]^2 + 0.75 * s[t - 1]^2 + as.numeric(v[t,] %*% xi))
    tail(e * s, n_obs)
}
y_vreg2 <- xts(cbind(sim_fun_mc(cbind(vr1, vr3), c(0.3, 0.08)), y_vreg[,2], y_vreg[,3]),
               order.by = index(y_vreg))
colnames(y_vreg2) <- colnames(y)[test_series]
vreg2_insample <- list(xts(cbind(vr1, vr3), order.by = index(y_vreg2)),
                       vreg_insample[[2]],
                       NULL)
x_vreg2 <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y_vreg2[,i], model = "garch", vreg = vreg2_insample[[i]])
    estimate(spec, keep_tmb = TRUE)
})
names(x_vreg2) <- colnames(y_vreg2)
x_vreg2 <- to_multi_estimate(x_vreg2)
dcc_vreg2_estimate <- estimate(dcc_modelspec(x_vreg2, dynamics = "dcc", distribution = "mvn"),
                               control = list(trace = 0), return_hessian = FALSE)

test_that("dcc predict propagates vreg to variance paths",{
    h <- 3
    nv_big <- list(matrix(10, nrow = h, ncol = 1), matrix(10, nrow = h, ncol = 1), NULL)
    names(nv_big) <- sn
    nv_zero <- list(matrix(0, nrow = h, ncol = 1), matrix(0, nrow = h, ncol = 1), NULL)
    names(nv_zero) <- sn
    p_big <- predict(dcc_vreg_estimate, h = h, nsim = 500, seed = 42, newvreg = nv_big)
    p_zero <- predict(dcc_vreg_estimate, h = h, nsim = 500, seed = 42, newvreg = nv_zero)
    H_big <- tscov(p_big, distribution = FALSE)
    H_zero <- tscov(p_zero, distribution = FALSE)
    # the variance paths differ materially for the series with vreg
    expect_true(max(abs(H_big[1,1,] - H_zero[1,1,])) > 0.01)
    expect_true(max(abs(H_big[2,2,] - H_zero[2,2,])) > 0.01)
    # and are identical for the series without vreg
    expect_equal(H_big[3,3,], H_zero[3,3,])
})

test_that("dcc predict variance level matches univariate analytic forecast",{
    # H[i,i,] is sigma_i^2 (R[i,i] == 1); the mean across draws of the simulated
    # variance equals tsgarch's analytic forecast, which pins the %*% xi level
    h <- 3
    p <- predict(dcc_vreg_estimate, h = h, nsim = 2000, seed = 42, newvreg = newvreg_hor)
    H <- tscov(p, distribution = FALSE)
    s1 <- as.numeric(predict(x_vreg[[1]], h = h, newvreg = newvreg_hor[[1]])$sigma)
    s2 <- as.numeric(predict(x_vreg[[2]], h = h, newvreg = newvreg_hor[[2]])$sigma)
    expect_equal(H[1,1,], s1^2, tolerance = 0.05)
    expect_equal(H[2,2,], s2^2, tolerance = 0.05)
})

test_that("dcc multi-column vreg propagates with correct ordering and level",{
    h <- 3
    nv_hor <- list(matrix(c(runif(h), 5 * runif(h)), ncol = 2), matrix(runif(h), ncol = 1), NULL)
    names(nv_hor) <- sn
    nv_zero <- list(matrix(0, nrow = h, ncol = 2), matrix(0, nrow = h, ncol = 1), NULL)
    names(nv_zero) <- sn
    p <- predict(dcc_vreg2_estimate, h = h, nsim = 2000, seed = 42, newvreg = nv_hor)
    p_zero <- predict(dcc_vreg2_estimate, h = h, nsim = 2000, seed = 42, newvreg = nv_zero)
    H <- tscov(p, distribution = FALSE)
    H_zero <- tscov(p_zero, distribution = FALSE)
    expect_true(max(abs(H[1,1,] - H_zero[1,1,])) > 0.01)
    expect_true(max(abs(H[2,2,] - H_zero[2,2,])) > 0.01)
    expect_equal(H[3,3,], H_zero[3,3,])
    # the level check against tsgarch's analytic forecast is what actually pins
    # the column ordering of the %*% xi pre-multiplication
    s1 <- as.numeric(predict(x_vreg2[[1]], h = h, newvreg = nv_hor[[1]])$sigma)
    s2 <- as.numeric(predict(x_vreg2[[2]], h = h, newvreg = nv_hor[[2]])$sigma)
    expect_equal(H[1,1,], s1^2, tolerance = 0.05)
    expect_equal(H[2,2,], s2^2, tolerance = 0.05)
    # simulate end to end with h + burn rows on the 2-column regressor
    nv_sim <- list(matrix(c(runif(5), 5 * runif(5)), ncol = 2), matrix(runif(5), ncol = 1), NULL)
    names(nv_sim) <- sn
    nv_sim_zero <- lapply(nv_sim, function(z) if (is.null(z)) NULL else z * 0)
    s_big <- simulate(dcc_vreg2_estimate, h = 3, burn = 2, nsim = 500, seed = 42, vreg = nv_sim)
    s_zero <- simulate(dcc_vreg2_estimate, h = 3, burn = 2, nsim = 500, seed = 42, vreg = nv_sim_zero)
    Hs <- tscov(s_big, distribution = FALSE)
    Hs_zero <- tscov(s_zero, distribution = FALSE)
    expect_true(max(abs(Hs[1,1,] - Hs_zero[1,1,])) > 0.01)
    expect_true(max(abs(Hs[2,2,] - Hs_zero[2,2,])) > 0.01)
    expect_equal(Hs[3,3,], Hs_zero[3,3,])
})

test_that("dcc simulate propagates vreg to variance paths",{
    nv_big <- list(matrix(10, nrow = 5, ncol = 1), matrix(10, nrow = 5, ncol = 1), NULL)
    names(nv_big) <- sn
    nv_zero <- list(matrix(0, nrow = 5, ncol = 1), matrix(0, nrow = 5, ncol = 1), NULL)
    names(nv_zero) <- sn
    s_big <- simulate(dcc_vreg_estimate, h = 3, burn = 2, nsim = 500, seed = 42, vreg = nv_big)
    s_zero <- simulate(dcc_vreg_estimate, h = 3, burn = 2, nsim = 500, seed = 42, vreg = nv_zero)
    H_big <- tscov(s_big, distribution = FALSE)
    H_zero <- tscov(s_zero, distribution = FALSE)
    expect_true(max(abs(H_big[1,1,] - H_zero[1,1,])) > 0.01)
    expect_true(max(abs(H_big[2,2,] - H_zero[2,2,])) > 0.01)
    expect_equal(H_big[3,3,], H_zero[3,3,])
})

test_that("dcc filter with named newvreg matches univariate filters",{
    filtered <- tsfilter(dcc_vreg_estimate, y = ynew, newvreg = newvreg_filter, update = FALSE)
    f_uni <- lapply(seq_along(test_series), function(i){
        tsfilter(x_vreg[[i]], y = ynew[,i], newvreg = newvreg_filter[[i]])
    })
    for (i in seq_along(test_series)) {
        expect_equal(coredata(sigma(filtered$spec$univariate[[i]])), coredata(sigma(f_uni[[i]])))
    }
    # materially different from zero regressors
    nv_zero <- newvreg_filter
    nv_zero[[1]] <- nv_zero[[1]] * 0
    nv_zero[[2]] <- nv_zero[[2]] * 0
    filtered_zero <- tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_zero, update = FALSE)
    expect_true(max(abs(coredata(sigma(filtered$spec$univariate)) - coredata(sigma(filtered_zero$spec$univariate)))) > 0.001)
})

test_that("missing newvreg is an error naming the affected series",{
    # whole argument NULL
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, update = FALSE), sn[1])
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, update = FALSE), sn[2])
    # NULL element for a series with vreg
    nv_null <- newvreg_filter
    nv_null[2] <- list(NULL)
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_null, update = FALSE), sn[2])
    # partial named list omitting a series with vreg
    nv_partial <- newvreg_filter[1]
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_partial, update = FALSE), sn[2])
    # same policy in predict and simulate
    expect_error(predict(dcc_vreg_estimate, h = 3, nsim = 10, seed = 1), sn[1])
    expect_error(simulate(dcc_vreg_estimate, h = 3, nsim = 10, seed = 1), sn[1])
})

test_that("vreg validation errors and warnings",{
    # bare matrix is not accepted
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = matrix(runif(10), ncol = 1), update = FALSE), "list")
    # unnamed list of the wrong length
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = list(matrix(runif(10), ncol = 1)), update = FALSE), "one element per series")
    # named list with unknown name
    nv_unknown <- newvreg_filter
    names(nv_unknown)[1] <- "UNKNOWN"
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_unknown, update = FALSE), "not matching")
    # wrong ncol
    nv_ncol <- newvreg_filter
    nv_ncol[[1]] <- matrix(runif(20), ncol = 2)
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_ncol, update = FALSE), "columns")
    # non-finite values
    nv_na <- newvreg_filter
    nv_na[[1]][1,1] <- NA
    expect_error(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_na, update = FALSE), sn[1])
    # element supplied for series without regressors warns and is ignored
    nv_extra <- newvreg_filter
    nv_extra[[3]] <- matrix(runif(10), ncol = 1)
    expect_warning(tsfilter(dcc_vreg_estimate, y = ynew, newvreg = nv_extra, update = FALSE), "Ignoring")
    # simulate requires h + burn rows
    nv_bad <- lapply(newvreg_sim, function(x) if (is.null(x)) NULL else x[1:3,,drop = FALSE])
    expect_error(simulate(dcc_vreg_estimate, h = 3, burn = 2, nsim = 50, seed = 42, vreg = nv_bad), "rows")
})

test_that("multiplicative vreg propagates end to end",{
    x_mult <- lapply(seq_along(test_series), function(i){
        spec <- garch_modelspec(y_vreg[,i], model = "garch", vreg = vreg_insample[[i]], multiplicative = TRUE)
        estimate(spec, keep_tmb = TRUE)
    })
    names(x_mult) <- colnames(y_vreg)
    x_mult <- to_multi_estimate(x_mult)
    mod <- estimate(dcc_modelspec(x_mult, dynamics = "dcc", distribution = "mvn"), control = list(trace = 0), return_hessian = FALSE)
    h <- 3
    nv_big <- list(matrix(10, nrow = h, ncol = 1), matrix(10, nrow = h, ncol = 1), NULL)
    names(nv_big) <- sn
    nv_zero <- list(matrix(0, nrow = h, ncol = 1), matrix(0, nrow = h, ncol = 1), NULL)
    names(nv_zero) <- sn
    p_big <- predict(mod, h = h, nsim = 500, seed = 42, newvreg = nv_big)
    p_zero <- predict(mod, h = h, nsim = 500, seed = 42, newvreg = nv_zero)
    H_big <- tscov(p_big, distribution = FALSE)
    H_zero <- tscov(p_zero, distribution = FALSE)
    expect_true(max(abs(H_big[1,1,] - H_zero[1,1,])) > 0.01)
    expect_true(max(abs(H_big[2,2,] - H_zero[2,2,])) > 0.01)
    expect_equal(H_big[3,3,], H_zero[3,3,])
})

test_that("cgarch predict propagates vreg to variance paths",{
    h <- 3
    nv_big <- list(matrix(10, nrow = h, ncol = 1), matrix(10, nrow = h, ncol = 1), NULL)
    names(nv_big) <- sn
    nv_zero <- list(matrix(0, nrow = h, ncol = 1), matrix(0, nrow = h, ncol = 1), NULL)
    names(nv_zero) <- sn
    p_big <- predict(cgarch_vreg_estimate, h = h, nsim = 500, seed = 42, newvreg = nv_big)
    p_zero <- predict(cgarch_vreg_estimate, h = h, nsim = 500, seed = 42, newvreg = nv_zero)
    H_big <- tscov(p_big, distribution = FALSE)
    H_zero <- tscov(p_zero, distribution = FALSE)
    expect_true(max(abs(H_big[1,1,] - H_zero[1,1,])) > 0.01)
    expect_true(max(abs(H_big[2,2,] - H_zero[2,2,])) > 0.01)
    expect_equal(H_big[3,3,], H_zero[3,3,])
})

test_that("cgarch simulate propagates vreg to variance paths",{
    nv_big <- list(matrix(10, nrow = 5, ncol = 1), matrix(10, nrow = 5, ncol = 1), NULL)
    names(nv_big) <- sn
    nv_zero <- list(matrix(0, nrow = 5, ncol = 1), matrix(0, nrow = 5, ncol = 1), NULL)
    names(nv_zero) <- sn
    s_big <- simulate(cgarch_vreg_estimate, h = 3, burn = 2, nsim = 500, seed = 42, vreg = nv_big)
    s_zero <- simulate(cgarch_vreg_estimate, h = 3, burn = 2, nsim = 500, seed = 42, vreg = nv_zero)
    H_big <- tscov(s_big, distribution = FALSE)
    H_zero <- tscov(s_zero, distribution = FALSE)
    expect_true(max(abs(H_big[1,1,] - H_zero[1,1,])) > 0.01)
    expect_true(max(abs(H_big[2,2,] - H_zero[2,2,])) > 0.01)
    expect_equal(H_big[3,3,], H_zero[3,3,])
})

test_that("cgarch filter with named newvreg matches univariate filters",{
    filtered <- tsfilter(cgarch_vreg_estimate, y = ynew, newvreg = newvreg_filter, update = FALSE)
    f_uni <- lapply(seq_along(test_series), function(i){
        tsfilter(x_vreg[[i]], y = ynew[,i], newvreg = newvreg_filter[[i]])
    })
    for (i in seq_along(test_series)) {
        expect_equal(coredata(sigma(filtered$spec$univariate[[i]])), coredata(sigma(f_uni[[i]])))
    }
})

test_that("cgarch missing newvreg is an error naming the affected series",{
    expect_error(tsfilter(cgarch_vreg_estimate, y = ynew, update = FALSE), sn[1])
    nv_partial <- newvreg_filter[2]
    expect_error(tsfilter(cgarch_vreg_estimate, y = ynew, newvreg = nv_partial, update = FALSE), sn[1])
    expect_error(predict(cgarch_vreg_estimate, h = 3, nsim = 10, seed = 1), sn[1])
    expect_error(simulate(cgarch_vreg_estimate, h = 3, nsim = 10, seed = 1), sn[1])
})
