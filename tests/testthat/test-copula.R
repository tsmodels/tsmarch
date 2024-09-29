test_that("cgarch constant residuals",{
    stdize <- FALSE
    type <- "standard"
    r <- coredata(residuals(global_cgarch_constant_estimate_e, standardize = stdize, type = type))
    expect_equal(as.numeric(r), as.numeric(global_cgarch_constant_estimate_e$spec$target$y))
    stdize <- TRUE
    r <- coredata(residuals(global_cgarch_constant_estimate_e, standardize = stdize, type = type))
    sd_r <- apply(r, 2, sd)
    expect_equal(mean(sd_r), 1.0, tolerance = 0.02)
})

test_that("cgarch dynamic residuals",{
    stdize <- FALSE
    type <- "standard"
    r <- coredata(residuals(global_cgarch_constant_estimate_e, standardize = stdize, type = type))
    expect_equal(as.numeric(r), as.numeric(global_cgarch_constant_estimate_e$spec$target$y))
    stdize <- TRUE
    r <- coredata(residuals(global_cgarch_constant_estimate_e, standardize = stdize, type = type))
    # these are the same as the constant case since standardization is on the univariate residuals (same model)
    sd_r <- apply(r, 2, sd)
    expect_equal(mean(sd_r), 1.0, tolerance = 0.02)
})

test_that("cgarch constant filter no update",{
    global_cgarch_constant_filter_p <- tsfilter(global_cgarch_constant_estimate_p, y = y[test_index,test_series], update = FALSE)
    global_cgarch_constant_filter_spd <- tsfilter(global_cgarch_constant_estimate_spd, y = y[test_index,test_series], update = FALSE)
    global_cgarch_constant_filter_e <- tsfilter(global_cgarch_constant_estimate_e, y = y[test_index,test_series], update = FALSE)

    expect_equal(residuals(global_cgarch_constant_estimate_p, type = "model"), residuals(global_cgarch_constant_filter_p, type = "model")[sample_index,], tolerance = 1e-5)
    expect_equal(residuals(global_cgarch_constant_estimate_spd, type = "model"), residuals(global_cgarch_constant_filter_spd, type = "model")[sample_index,], tolerance = 1e-5)
    expect_equal(residuals(global_cgarch_constant_estimate_e, type = "model"), residuals(global_cgarch_constant_filter_e, type = "model")[sample_index,], tolerance = 1e-5)
})

test_that("cgarch dynamic filter no update",{
    global_cgarch_dcc_filter_p <- tsfilter(global_cgarch_dcc_estimate_p, y = y[test_index,test_series], update = FALSE)
    global_cgarch_dcc_filter_spd <- tsfilter(global_cgarch_dcc_estimate_spd, y = y[test_index,test_series], update = FALSE)
    global_cgarch_dcc_filter_e <- tsfilter(global_cgarch_dcc_estimate_e, y = y[test_index,test_series], update = FALSE)

    expect_equal(residuals(global_cgarch_dcc_estimate_p, type = "model"), residuals(global_cgarch_dcc_filter_p, type = "model")[sample_index,], tolerance = 1e-8)
    expect_equal(residuals(global_cgarch_dcc_estimate_spd, type = "model"), residuals(global_cgarch_dcc_filter_spd, type = "model")[sample_index,], tolerance = 1e-8)
    expect_equal(residuals(global_cgarch_dcc_estimate_e, type = "model"), residuals(global_cgarch_dcc_filter_e, type = "model")[sample_index,], tolerance = 1e-8)
})

test_that("cgarch dynamic prediction",{
    n_series <- length(test_series)
    h <- 1
    nsim <- 10
    p <- predict(global_cgarch_dcc_estimate_p, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(n_series, n_series, h, nsim))
    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))

    h <- 3
    p <- predict(global_cgarch_dcc_estimate_p, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(n_series, n_series, h, nsim))
    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))

})

test_that("cgarch constant prediction",{
    n_series <- length(test_series)
    h <- 1
    nsim <- 10
    p <- predict(global_cgarch_constant_estimate_spd, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    # constant correlation no distribution
    expect_equal(dim(C), c(n_series, n_series))
    C <- tscor(p, distribution = FALSE)
    # constant correlation no distribution
    expect_equal(dim(C), c(n_series, n_series))

    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))

    h <- 3
    p <- predict(global_cgarch_constant_estimate_spd, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(n_series, n_series))
    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))
})
