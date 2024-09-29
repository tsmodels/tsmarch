test_that("dcc constant residuals",{
    stdize <- FALSE
    type <- "standard"
    r <- coredata(residuals(global_dcc_constant_estimate, standardize = stdize, type = type))
    expect_equal(as.numeric(r), as.numeric(global_dcc_constant_estimate$spec$target$y))
    stdize <- TRUE
    r <- coredata(residuals(global_dcc_constant_estimate, standardize = stdize, type = type))
    sd_r <- apply(r, 2, sd)
    expect_equal(mean(sd_r), 1.0, tolerance = 0.02)
})

test_that("dcc dynamic residuals",{
    stdize <- FALSE
    type <- "standard"
    r <- coredata(residuals(global_dcc_dynamic_estimate, standardize = stdize, type = type))
    expect_equal(as.numeric(r), as.numeric(global_dcc_dynamic_estimate$spec$target$y))
    stdize <- TRUE
    r <- coredata(residuals(global_dcc_dynamic_estimate, standardize = stdize, type = type))
    # these are the same as the constant case since standardization is on the univariate residuals (same model)
    sd_r <- apply(r, 2, sd)
    expect_equal(mean(sd_r), 1.0, tolerance = 0.02)
})


test_that("dcc constant filter no update",{
    global_dcc_constant_filter <- tsfilter(global_dcc_constant_estimate, y = y[test_index,test_series], update = FALSE)
    expect_equal(residuals(global_dcc_constant_estimate, standardize = TRUE), residuals(global_dcc_constant_filter, standardize = TRUE)[sample_index,], tolerance = 1e-5)
})

test_that("dcc dynamic filter no update",{
    global_dcc_dynamic_filter <- tsfilter(global_dcc_dynamic_estimate, y = y[test_index,test_series], update = FALSE)
    global_adcc_dynamic_filter <- tsfilter(global_adcc_dynamic_estimate, y = y[test_index,test_series], update = FALSE)
    expect_equal(residuals(global_dcc_dynamic_estimate, standardize = TRUE), residuals(global_dcc_dynamic_filter, standardize = TRUE)[sample_index,], tolerance = 1e-5)
    expect_equal(residuals(global_adcc_dynamic_estimate, standardize = TRUE), residuals(global_adcc_dynamic_filter, standardize = TRUE)[sample_index,], tolerance = 1e-5)
})

test_that("dcc dynamic prediction",{
    n_series <- length(test_series)
    h <- 1
    nsim <- 10
    p <- predict(global_dcc_dynamic_estimate, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(n_series, n_series, h, nsim))
    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))

    h <- 3
    p <- predict(global_dcc_dynamic_estimate, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(n_series, n_series, h, nsim))
    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))
})

test_that("dcc constant prediction",{
    n_series <- length(test_series)
    h <- 1
    nsim <- 10
    p <- predict(global_dcc_constant_estimate, h = h, nsim = nsim, seed = 100)
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
    p <- predict(global_dcc_constant_estimate, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(n_series, n_series, h, nsim))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(n_series, n_series))
    port <- tsaggregate(p, weights = rep(1/n_series, n_series))
    expect_equal(dim(port$mu), c(nsim, h))
    expect_equal(dim(port$sigma), c(nsim, h))
})
