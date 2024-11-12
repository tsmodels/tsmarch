test_that("GOGARCH Estimation Convolution/Geometric Moments",{
    series <- 6
    geometric_pmoments <- tsaggregate(global_gogarch_mod_dr, weights = rep(1/series, series))
    cf <- tsconvolve(global_gogarch_mod_dr, weights = rep(1/series, series), fft_support = c(-1, 1))
    n <- length(gogarch_sample_index)
    int_moments <- matrix(0, ncol = 2, nrow = n)
    for (i in 1:n) {
        # we've already fixed the support to (-1, 1), but show here that we can also
        # retrieve the support range if it instead calculated
        rnge <- attr(cf$y[[i]], "support")
        dx <- dfft(cf, index = i)
        f1 <- function(x) x * dx(x)
        f_mu <- integrate(f1, rnge[1], rnge[2])$value
        f3 <- function(x) (x - f_mu)^3 * dx(x)
        f4 <- function(x) (x - f_mu)^4 * dx(x)
        fx3 <- integrate(f3, rnge[1], rnge[2])$value
        fx4 <- integrate(f4, rnge[1], rnge[2])$value
        int_moments[i,1] <- fx3
        int_moments[i,2] <- fx4
    }
    int_moments[,1] <- int_moments[,1]/as.numeric(geometric_pmoments$sigma)^3
    int_moments[,2] <- int_moments[,2]/as.numeric(geometric_pmoments$sigma)^4
    expect_equal(as.numeric(geometric_pmoments$skewness), int_moments[,1], tolerance = 0.1)
    expect_equal(as.numeric(geometric_pmoments$kurtosis), int_moments[,2], tolerance = 0.1)
})

test_that("GOGARCH Prediction Convolution/Geometric Moments [h=1]",{
    series <- 6
    nsim <- 10
    h <- 1
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    geometric_pmoments <- tsaggregate(p, weights = rep(1/series, series))
    cf <- tsconvolve(p, weights = rep(1/series, series), fft_support = c(-1, 1), distribution = TRUE)
    kurtosis_moment <- skew_moment <- matrix(0, ncol = h, nrow = nsim)
    for (i in 1:h) {
        for (j in 1:nsim) {
            rnge <- attr(cf$y[[j]][[i]], "support")
            dx <- dfft(cf, index = i, sim = j)
            f1 <- function(x) x * dx(x)
            f_mu <- integrate(f1, rnge[1], rnge[2])$value
            f3 <- function(x) (x - f_mu)^3 * dx(x)
            f4 <- function(x) (x - f_mu)^4 * dx(x)
            fx3 <- integrate(f3, rnge[1], rnge[2])$value
            fx4 <- integrate(f4, rnge[1], rnge[2])$value
            skew_moment[j,i] <- fx3/as.numeric(geometric_pmoments$sigma[j,i])^3
            kurtosis_moment[j,i] <- fx4/as.numeric(geometric_pmoments$sigma[j,i])^4
        }
    }
    expect_equal(as.numeric(skew_moment), as.numeric(geometric_pmoments$skewness), tolerance = 0.1)
    expect_equal(as.numeric(kurtosis_moment), as.numeric(geometric_pmoments$kurtosis), tolerance = 0.1)
})

test_that("GOGARCH Prediction Convolution/Geometric Moments [h>1]",{

    series <- 6
    nsim <- 10
    h <- 3
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    geometric_pmoments <- tsaggregate(p, weights = rep(1/series, series))
    cf <- tsconvolve(p, weights = rep(1/series, series), fft_support = c(-1, 1), distribution = TRUE)
    kurtosis_moment <- skew_moment <- matrix(0, ncol = h, nrow = nsim)
    for (i in 1:h) {
        for (j in 1:nsim) {
            rnge <- attr(cf$y[[j]][[i]], "support")
            dx <- dfft(cf, index = i, sim = j)
            f1 <- function(x) x * dx(x)
            f_mu <- integrate(f1, rnge[1], rnge[2])$value
            f3 <- function(x) (x - f_mu)^3 * dx(x)
            f4 <- function(x) (x - f_mu)^4 * dx(x)
            fx3 <- integrate(f3, rnge[1], rnge[2])$value
            fx4 <- integrate(f4, rnge[1], rnge[2])$value
            skew_moment[j,i] <- fx3/as.numeric(geometric_pmoments$sigma[j,i])^3
            kurtosis_moment[j,i] <- fx4/as.numeric(geometric_pmoments$sigma[j,i])^4
        }
    }
    expect_equal(as.numeric(skew_moment), as.numeric(geometric_pmoments$skewness), tolerance = 0.1)
    expect_equal(as.numeric(kurtosis_moment), as.numeric(geometric_pmoments$kurtosis), tolerance = 0.1)
})


test_that("GOGARCH Simulation Convolution/Geometric Moments [h=1]",{
    series <- 6
    nsim <- 10
    h <- 1
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    geometric_pmoments <- tsaggregate(p, weights = rep(1/series, series))
    cf <- tsconvolve(p, weights = rep(1/series, series), fft_support = c(-1, 1), distribution = TRUE)
    kurtosis_moment <- skew_moment <- matrix(0, ncol = h, nrow = nsim)
    for (i in 1:h) {
        for (j in 1:nsim) {
            rnge <- attr(cf$y[[j]][[i]], "support")
            dx <- dfft(cf, index = i, sim = j)
            f1 <- function(x) x * dx(x)
            f_mu <- integrate(f1, rnge[1], rnge[2])$value
            f3 <- function(x) (x - f_mu)^3 * dx(x)
            f4 <- function(x) (x - f_mu)^4 * dx(x)
            fx3 <- integrate(f3, rnge[1], rnge[2])$value
            fx4 <- integrate(f4, rnge[1], rnge[2])$value
            skew_moment[j,i] <- fx3/as.numeric(geometric_pmoments$sigma[j,i])^3
            kurtosis_moment[j,i] <- fx4/as.numeric(geometric_pmoments$sigma[j,i])^4
        }
    }
    expect_equal(as.numeric(skew_moment), as.numeric(geometric_pmoments$skewness), tolerance = 0.1)
    expect_equal(as.numeric(kurtosis_moment), as.numeric(geometric_pmoments$kurtosis), tolerance = 0.1)
})

test_that("GOGARCH Simulation Convolution/Geometric Moments [h>1]",{
    series <- 6
    nsim <- 10
    h <- 3
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    geometric_pmoments <- tsaggregate(p, weights = rep(1/series, series))
    cf <- tsconvolve(p, weights = rep(1/series, series), fft_support = c(-1, 1), distribution = TRUE)
    kurtosis_moment <- skew_moment <- matrix(0, ncol = h, nrow = nsim)
    for (i in 1:h) {
        for (j in 1:nsim) {
            rnge <- attr(cf$y[[j]][[i]], "support")
            dx <- dfft(cf, index = i, sim = j)
            f1 <- function(x) x * dx(x)
            f_mu <- integrate(f1, rnge[1], rnge[2], abs.tol = 1e-8, stop.on.error = FALSE)$value
            f3 <- function(x) (x - f_mu)^3 * dx(x)
            f4 <- function(x) (x - f_mu)^4 * dx(x)
            fx3 <- integrate(f3, rnge[1], rnge[2])$value
            fx4 <- integrate(f4, rnge[1], rnge[2])$value
            skew_moment[j,i] <- fx3/as.numeric(geometric_pmoments$sigma[j,i])^3
            kurtosis_moment[j,i] <- fx4/as.numeric(geometric_pmoments$sigma[j,i])^4
        }
    }
    expect_equal(as.numeric(skew_moment), as.numeric(geometric_pmoments$skewness), tolerance = 0.1)
    expect_equal(as.numeric(kurtosis_moment), as.numeric(geometric_pmoments$kurtosis), tolerance = 0.1)
})
