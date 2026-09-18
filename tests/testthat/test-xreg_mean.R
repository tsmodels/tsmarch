# fixtures: small and fast, built locally (helper-global.R keeps the non-xreg globals)
# series with genuine regressor effects in the mean so the effect is material
set.seed(1042)
n_obs <- 500
xr1 <- rnorm(n_obs)
xr2 <- 0.5 * xr1 + rnorm(n_obs, sd = 0.8)
xr3 <- rnorm(n_obs)
sim_fun <- function(xreg, tau) {
    e <- rnorm(n_obs + 50)
    s <- rep(1, n_obs + 50)
    for (t in 2:(n_obs + 50)) s[t] <- sqrt(0.05 + 0.1 * e[t - 1]^2 + 0.85 * s[t - 1]^2)
    z <- e * s
    if (is.null(xreg)) return(tail(z, n_obs))
    tail(tau * c(rep(0, 50), xreg) + z, n_obs)
}
y_xreg <- xts(cbind(sim_fun(xr1, 0.8), sim_fun(xr2, -0.6), sim_fun(NULL, 0)),
              order.by = y[1:n_obs, test_series] |> index())
colnames(y_xreg) <- colnames(y)[test_series]
xreg_insample <- list(xts(xr1, order.by = index(y_xreg)),
                      xts(xr2, order.by = index(y_xreg)),
                      NULL)

x_xonly <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y_xreg[,i], model = "garch", arma = c(0,0), xreg = xreg_insample[[i]])
    estimate(spec, keep_tmb = TRUE)
})
names(x_xonly) <- colnames(y_xreg)
x_xonly <- to_multi_estimate(x_xonly)

x_armax <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y_xreg[,i], model = "garch", arma = c(1,0), xreg = xts(get(paste0("xr", i)), order.by = index(y_xreg)))
    estimate(spec, keep_tmb = TRUE)
})
names(x_armax) <- colnames(y_xreg)
x_armax <- to_multi_estimate(x_armax)

dcc_xreg_spec <- dcc_modelspec(x_xonly, dynamics = "dcc", distribution = "mvn")
dcc_xreg_estimate <- estimate(dcc_xreg_spec, control = list(trace = 0), return_hessian = FALSE)
cgarch_xreg_spec <- cgarch_modelspec(x_xonly, dynamics = "constant", transformation = "parametric", copula = "mvn")
cgarch_xreg_estimate <- estimate(cgarch_xreg_spec, control = list(trace = 0))

sn <- colnames(y_xreg)
new_dates <- seq(tail(index(y), 1), length.out = 10, by = 7)
ynew <- xts(matrix(rnorm(30), ncol = 3), order.by = new_dates)
colnames(ynew) <- sn
newxreg_filter <- list(matrix(rnorm(10), ncol = 1), matrix(rnorm(10), ncol = 1), NULL)
names(newxreg_filter) <- sn
newxreg_hor <- list(matrix(rnorm(3), ncol = 1), matrix(rnorm(3), ncol = 1), NULL)
names(newxreg_hor) <- sn
newxreg_sim <- list(matrix(rnorm(5), ncol = 1), matrix(rnorm(5), ncol = 1), NULL)
names(newxreg_sim) <- sn

test_that("xreg mean propagates to dcc spec",{
    expect_equal(dcc_xreg_spec$target$mu, coredata(fitted(x_xonly)))
    per_column_range <- apply(dcc_xreg_spec$target$mu, 2, function(z) diff(range(z)))
    # the regressor effect is present in-sample for the series with regressors
    expect_true(all(per_column_range[1:2] > 0))
})

test_that("xreg mean propagates to cgarch spec",{
    expect_equal(cgarch_xreg_spec$target$mu, coredata(fitted(x_xonly)))
    per_column_range <- apply(cgarch_xreg_spec$target$mu, 2, function(z) diff(range(z)))
    expect_true(all(per_column_range[1:2] > 0))
})

test_that("cond_mean rejected at specification with xreg first stage",{
    m_spec <- matrix(0, ncol = length(test_series), nrow = NROW(y_xreg))
    expect_error(dcc_modelspec(x_xonly, cond_mean = m_spec), "regressor")
    expect_error(cgarch_modelspec(x_xonly, cond_mean = m_spec), "regressor")
})

test_that("cond_mean rejected in filter simulate and predict",{
    m_new <- matrix(0, ncol = length(test_series), nrow = NROW(ynew))
    m_hor <- matrix(0, ncol = length(test_series), nrow = 3)
    expect_error(tsfilter(dcc_xreg_estimate, y = ynew, cond_mean = m_new), "regressor")
    expect_error(simulate(dcc_xreg_estimate, h = 3, nsim = 5, cond_mean = m_hor), "regressor")
    expect_error(predict(dcc_xreg_estimate, h = 3, nsim = 5, cond_mean = m_hor), "regressor")
    expect_error(tsfilter(cgarch_xreg_estimate, y = ynew, cond_mean = m_new), "regressor")
    expect_error(simulate(cgarch_xreg_estimate, h = 3, nsim = 5, cond_mean = m_hor), "regressor")
    expect_error(predict(cgarch_xreg_estimate, h = 3, nsim = 5, cond_mean = m_hor), "regressor")
})

test_that("dcc filter with named newxreg matches univariate filters",{
    filtered <- tsfilter(dcc_xreg_estimate, y = ynew, newxreg = newxreg_filter, update = FALSE)
    expect_equal(coredata(fitted(filtered)), coredata(fitted(filtered$spec$univariate)))
    f_uni <- lapply(seq_along(test_series), function(i){
        tsfilter(x_xonly[[i]], y = ynew[,i], newxreg = newxreg_filter[[i]])
    })
    expect_equal(unname(coredata(fitted(filtered))), unname(do.call(cbind, lapply(f_uni, function(f) coredata(fitted(f))))))
    # materially different from zero regressors
    nx_zero <- newxreg_filter
    nx_zero[[1]] <- nx_zero[[1]] * 0
    nx_zero[[2]] <- nx_zero[[2]] * 0
    filtered_zero <- tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_zero, update = FALSE)
    expect_true(max(abs(coredata(fitted(filtered)) - coredata(fitted(filtered_zero)))) > 0.01)
})

test_that("cgarch filter with named newxreg matches univariate filters",{
    filtered <- tsfilter(cgarch_xreg_estimate, y = ynew, newxreg = newxreg_filter, update = FALSE)
    expect_equal(coredata(fitted(filtered)), coredata(fitted(filtered$spec$univariate)))
    f_uni <- lapply(seq_along(test_series), function(i){
        tsfilter(x_xonly[[i]], y = ynew[,i], newxreg = newxreg_filter[[i]])
    })
    expect_equal(unname(coredata(fitted(filtered))), unname(do.call(cbind, lapply(f_uni, function(f) coredata(fitted(f))))))
})

test_that("dcc predict propagates xreg conditional mean",{
    h <- 3
    nsim <- 2000
    p <- predict(dcc_xreg_estimate, h = h, nsim = nsim, seed = 42, newxreg = newxreg_hor)
    mu_mean <- apply(p$mu, c(1,2), mean)
    mu_analytic <- do.call(cbind, lapply(seq_along(test_series), function(i){
        as.numeric(predict(x_xonly[[i]], h = h, newxreg = newxreg_hor[[i]])$mean)
    }))
    expect_equal(unname(mu_mean), unname(mu_analytic), tolerance = 0.1)
    # shifting the regressor by a constant shifts the mean by approximately tau * delta
    delta <- 1
    nx_shift <- newxreg_hor
    nx_shift[[1]] <- nx_shift[[1]] + delta
    p_shift <- predict(dcc_xreg_estimate, h = h, nsim = nsim, seed = 42, newxreg = nx_shift)
    mu_shift <- apply(p_shift$mu, c(1,2), mean)
    tau1 <- x_xonly[[1]]$parmatrix[parameter == "tau1"]$value
    expect_equal(unname(mu_shift[,1] - mu_mean[,1]), rep(tau1 * delta, h), tolerance = 0.1)
})

test_that("cgarch predict propagates xreg conditional mean",{
    h <- 3
    nsim <- 2000
    p <- predict(cgarch_xreg_estimate, h = h, nsim = nsim, seed = 42, newxreg = newxreg_hor)
    mu_mean <- apply(p$mu, c(1,2), mean)
    mu_analytic <- do.call(cbind, lapply(seq_along(test_series), function(i){
        as.numeric(predict(x_xonly[[i]], h = h, newxreg = newxreg_hor[[i]])$mean)
    }))
    expect_equal(unname(mu_mean), unname(mu_analytic), tolerance = 0.1)
    delta <- 1
    nx_shift <- newxreg_hor
    nx_shift[[1]] <- nx_shift[[1]] + delta
    p_shift <- predict(cgarch_xreg_estimate, h = h, nsim = nsim, seed = 42, newxreg = nx_shift)
    mu_shift <- apply(p_shift$mu, c(1,2), mean)
    tau1 <- x_xonly[[1]]$parmatrix[parameter == "tau1"]$value
    expect_equal(unname(mu_shift[,1] - mu_mean[,1]), rep(tau1 * delta, h), tolerance = 0.1)
})

test_that("dcc simulate with xreg requires h + burn rows",{
    nx_bad <- lapply(newxreg_sim, function(x) if (is.null(x)) NULL else x[1:3,,drop = FALSE])
    expect_error(simulate(dcc_xreg_estimate, h = 3, burn = 2, nsim = 500, seed = 42, xreg = nx_bad), "rows")
    s <- simulate(dcc_xreg_estimate, h = 3, burn = 2, nsim = 500, seed = 42, xreg = newxreg_sim)
    expect_equal(dim(s$mu), c(3, length(test_series), 500))
})

test_that("xreg validation errors and warnings",{
    # bare matrix is not accepted
    expect_error(tsfilter(dcc_xreg_estimate, y = ynew, newxreg = matrix(rnorm(10), ncol = 1), update = FALSE), "list")
    # unnamed list of the wrong length
    expect_error(tsfilter(dcc_xreg_estimate, y = ynew, newxreg = list(matrix(rnorm(10), ncol = 1)), update = FALSE), "one element per series")
    # named list with unknown name
    nx_unknown <- newxreg_filter
    names(nx_unknown)[1] <- "UNKNOWN"
    expect_error(tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_unknown, update = FALSE), "not matching")
    # wrong ncol
    nx_ncol <- newxreg_filter
    nx_ncol[[1]] <- matrix(rnorm(20), ncol = 2)
    expect_error(tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_ncol, update = FALSE), "columns")
    # non-finite values
    nx_na <- newxreg_filter
    nx_na[[1]][1,1] <- NA
    expect_error(tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_na, update = FALSE), sn[1])
    # element supplied for series without regressors warns and is ignored
    nx_extra <- newxreg_filter
    nx_extra[[3]] <- matrix(rnorm(10), ncol = 1)
    expect_warning(tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_extra, update = FALSE), "Ignoring")
    # partial named list omitting a series with regressors: single warning, zero-filled
    nx_partial <- newxreg_filter[1]
    expect_warning(
        filtered_partial <- tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_partial, update = FALSE),
        "Setting to zero")
    nx_explicit <- newxreg_filter
    nx_explicit[[2]] <- matrix(0, nrow = NROW(ynew), ncol = 1)
    filtered_explicit <- tsfilter(dcc_xreg_estimate, y = ynew, newxreg = nx_explicit, update = FALSE)
    expect_equal(coredata(fitted(filtered_partial)), coredata(fitted(filtered_explicit)))
})

test_that("armax first stage works end to end",{
    spec <- dcc_modelspec(x_armax, dynamics = "dcc", distribution = "mvn")
    mod <- estimate(spec, control = list(trace = 0), return_hessian = FALSE)
    expect_s3_class(mod, "dcc.estimate")
    nx <- lapply(seq_along(test_series), function(i) matrix(rnorm(NROW(ynew)), ncol = 1))
    filtered <- tsfilter(mod, y = ynew, newxreg = nx, update = FALSE)
    f_uni <- lapply(seq_along(test_series), function(i){
        tsfilter(x_armax[[i]], y = ynew[,i], newxreg = nx[[i]])
    })
    expect_equal(unname(coredata(fitted(filtered))), unname(do.call(cbind, lapply(f_uni, function(f) coredata(fitted(f))))))
    nx_hor <- lapply(seq_along(test_series), function(i) matrix(rnorm(3), ncol = 1))
    p <- predict(mod, h = 3, nsim = 50, seed = 42, newxreg = nx_hor)
    expect_equal(dim(p$mu), c(3, length(test_series), 50))
    nx_sim <- lapply(seq_along(test_series), function(i) matrix(rnorm(5), ncol = 1))
    s <- simulate(mod, h = 3, burn = 2, nsim = 50, seed = 42, xreg = nx_sim)
    expect_equal(dim(s$mu), c(3, length(test_series), 50))
})
