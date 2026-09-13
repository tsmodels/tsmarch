# fixtures: small and fast, built locally (helper-global.R keeps the non-ARMA globals)
# series with genuine AR(1) mean dynamics so the conditional mean is material
set.seed(1042)
n_obs <- 500
sim_fun <- function(phi) {
    e <- rnorm(n_obs + 50)
    s <- rep(1, n_obs + 50)
    for (t in 2:(n_obs + 50)) s[t] <- sqrt(0.05 + 0.1 * e[t - 1]^2 + 0.85 * s[t - 1]^2)
    z <- e * s
    yv <- numeric(n_obs + 50)
    for (t in 2:(n_obs + 50)) yv[t] <- phi * yv[t - 1] + z[t]
    tail(yv, n_obs)
}
y_arma <- xts(do.call(cbind, lapply(c(0.6, 0.4, 0.5), sim_fun)), order.by = y[1:n_obs, test_series] |> index())
colnames(y_arma) <- colnames(y)[test_series]
new_index <- tail(y[, test_series], 5) |> index()

x_arma <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y_arma[,i], model = "garch", arma = c(1,0))
    estimate(spec, keep_tmb = TRUE)
})
names(x_arma) <- colnames(y_arma)
x_arma <- to_multi_estimate(x_arma)

x_mixed <- lapply(seq_along(test_series), function(i){
    arma <- if (i == 1) c(1,0) else c(0,0)
    spec <- garch_modelspec(y_arma[,i], model = "garch", arma = arma)
    estimate(spec, keep_tmb = TRUE)
})
names(x_mixed) <- colnames(y_arma)
x_mixed <- to_multi_estimate(x_mixed)

x_arma3 <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y_arma[,i], model = "garch", order = c(1,1), arma = c(3,0))
    estimate(spec, keep_tmb = TRUE)
})
names(x_arma3) <- colnames(y_arma)
x_arma3 <- to_multi_estimate(x_arma3)

dcc_arma_spec <- dcc_modelspec(x_arma, dynamics = "dcc", distribution = "mvn")
dcc_arma_estimate <- estimate(dcc_arma_spec, control = list(trace = 0), return_hessian = FALSE)
cgarch_arma_spec <- cgarch_modelspec(x_arma, dynamics = "constant", transformation = "parametric", copula = "mvn")
cgarch_arma_estimate <- estimate(cgarch_arma_spec, control = list(trace = 0))

m_spec <- matrix(0, ncol = length(test_series), nrow = NROW(y_arma))
m_new <- matrix(0, ncol = length(test_series), nrow = 5)
m_hor <- matrix(0, ncol = length(test_series), nrow = 3)

test_that("arma mean propagates to dcc spec",{
    expect_equal(dcc_arma_spec$target$mu, coredata(fitted(x_arma)))
    per_column_range <- apply(dcc_arma_spec$target$mu, 2, function(z) diff(range(z)))
    expect_true(all(per_column_range > 0))
})

test_that("arma mean propagates to cgarch spec",{
    expect_equal(cgarch_arma_spec$target$mu, coredata(fitted(x_arma)))
    per_column_range <- apply(cgarch_arma_spec$target$mu, 2, function(z) diff(range(z)))
    expect_true(all(per_column_range > 0))
})

test_that("cond_mean rejected at specification with arma first stage",{
    expect_error(dcc_modelspec(x_arma, cond_mean = m_spec), "ARMA")
    expect_error(cgarch_modelspec(x_arma, cond_mean = m_spec), "ARMA")
    expect_error(dcc_modelspec(x_mixed, cond_mean = m_spec), "ARMA")
    expect_error(cgarch_modelspec(x_mixed, cond_mean = m_spec), "ARMA")
})

test_that("cond_mean rejected in filter simulate and predict",{
    ynew <- xts(matrix(rnorm(15), ncol = 3), order.by = new_index)
    colnames(ynew) <- colnames(y_arma)
    expect_error(tsfilter(dcc_arma_estimate, y = ynew, cond_mean = m_new), "ARMA")
    expect_error(simulate(dcc_arma_estimate, h = 3, nsim = 5, cond_mean = m_hor), "ARMA")
    expect_error(predict(dcc_arma_estimate, h = 3, nsim = 5, cond_mean = m_hor), "ARMA")
    expect_error(tsfilter(cgarch_arma_estimate, y = ynew, cond_mean = m_new), "ARMA")
    expect_error(simulate(cgarch_arma_estimate, h = 3, nsim = 5, cond_mean = m_hor), "ARMA")
    expect_error(predict(cgarch_arma_estimate, h = 3, nsim = 5, cond_mean = m_hor), "ARMA")
})

test_that("filter continues the arma recursion",{
    ynew <- xts(matrix(rnorm(30), ncol = 3), order.by = seq(tail(new_index, 1), length.out = 10, by = 7))
    colnames(ynew) <- colnames(y_arma)
    filtered <- tsfilter(dcc_arma_estimate, y = ynew, update = FALSE)
    expect_equal(coredata(fitted(filtered)), coredata(fitted(filtered$spec$univariate)))
    expect_equal(NROW(fitted(filtered)), NROW(y_arma) + NROW(ynew))
})

test_that("predict propagates arma conditional mean",{
    h <- 3
    nsim <- 2000
    p_arma <- predict(dcc_arma_estimate, h = h, nsim = nsim, seed = 42)
    mu_mean <- apply(p_arma$mu, c(1,2), mean)
    mu_analytic <- do.call(cbind, lapply(1:length(test_series), function(i) as.numeric(predict(x_arma[[i]], h = h)$mean)))
    # cross path mean matches the analytic arma point forecast
    expect_equal(unname(mu_mean), unname(mu_analytic), tolerance = 0.1)
    # and is materially different from zero (the arma mean is present)
    expect_true(max(abs(mu_mean)) > 0.1)
})

test_that("simulate propagates arma conditional mean",{
    # init_method = "end" continues the arma recursion from the last observations
    s <- simulate(dcc_arma_estimate, h = 3, nsim = 500, seed = 42, init_method = "end")
    expect_equal(dim(s$mu), c(3, length(test_series), 500))
    mu_mean <- apply(s$mu, c(1,2), mean)
    expect_true(max(abs(mu_mean)) > 0.1)
})

test_that("arma order greater than garch order works end to end",{
    spec <- dcc_modelspec(x_arma3, dynamics = "dcc", distribution = "mvn")
    mod <- estimate(spec, control = list(trace = 0), return_hessian = TRUE)
    expect_s3_class(mod, "dcc.estimate")
    ynew <- xts(matrix(rnorm(15), ncol = 3), order.by = new_index)
    colnames(ynew) <- colnames(y_arma)
    filtered <- tsfilter(mod, y = ynew, update = FALSE)
    expect_equal(coredata(fitted(filtered)), coredata(fitted(filtered$spec$univariate)))
    s <- simulate(mod, h = 2, nsim = 10, seed = 42)
    expect_equal(dim(s$mu), c(2, length(test_series), 10))
    p <- predict(mod, h = 2, nsim = 10, seed = 42)
    expect_equal(dim(p$mu), c(2, length(test_series), 10))
})

test_that("cond_mean still works without arma dynamics",{
    mu_in <- matrix(0.5, ncol = length(test_series), nrow = 1100)
    spec <- dcc_modelspec(x, cond_mean = mu_in)
    expect_equal(unname(spec$target$mu), unname(mu_in))
    mod <- estimate(spec, control = list(trace = 0), return_hessian = FALSE)
    expect_equal(unname(coredata(fitted(mod))), unname(mu_in))
    mu_new <- matrix(0.5, ncol = length(test_series), nrow = length(test_index))
    filtered <- tsfilter(mod, y = y[test_index, test_series], cond_mean = mu_new, update = FALSE)
    expect_equal(unname(coredata(fitted(filtered))[test_index,]), unname(mu_new))
    mu_hor <- matrix(0.5, ncol = length(test_series), nrow = 3)
    p <- predict(mod, h = 3, nsim = 100, cond_mean = mu_hor, seed = 42)
    expect_equal(unname(apply(p$mu, c(1,2), mean)), unname(mu_hor), tolerance = 0.2)
})
