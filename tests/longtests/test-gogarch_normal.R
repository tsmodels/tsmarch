test_that("GOGARCH Normal",{
    spec <- gogarch_modelspec(globalindices[gogarch_sample_index,1:n_series], distribution = "norm",
                              components = n_series)
    mod <- suppressWarnings(estimate(spec))
    a <- tsaggregate(mod, weights = rep(1/n_series, n_series))
    C <- tscov(mod)
    R <- tscor(mod)
    p <- predict(mod, h = 20, nsim = 2000, seed = 100)
    ap <- tsaggregate(p, weights = rep(1/n_series, n_series))
    Cp <- tscov(p)
    Rp <- tscor(p)

})
