# test fitted, tsaggregate, dfft, pfft and qfft

test_that("GOGARCH Estimation [cond_mean]",{
    w <- rep(1/6, 6)
    f <- coredata(fitted(mod_full))
    expect_equal(as.numeric(f), as.numeric(coredata(mu)))
    w_mu <- coredata(mu) %*% w
    a <- tsaggregate(mod_full, weights = rep(1/6, 6))
    expect_equal(as.numeric(w_mu), as.numeric(a$mu))
    cf <- tsconvolve(mod_full, weights = w, fft_support = c(-1, 1))
    # test dfft
    n <- length(gogarch_sample_index)
    test_mat <- matrix(0, ncol = 4, nrow = n)
    for (i in 1:n) {
        dx <- dfft(cf, index = i)
        f <- function(x) x * dx(x)
        w_imu <- integrate(f, -1, 1)$value
        test_mat[i,1] <- w_imu
        test_mat[i,2] <- w_mu[i]
        f3 <- function(x) (x - w_imu)^3 * dx(x)
        f4 <- function(x) (x - w_imu)^4 * dx(x)
        fx3 <- integrate(f3, -1, 1, abs.tol = 1e-8)$value
        fx4 <- integrate(f4, -1, 1, abs.tol = 1e-8)$value
        test_mat[i,3] <- fx3/as.numeric(a$sigma[i])^3
        test_mat[i,4] <- fx4/as.numeric(a$sigma[i])^4
    }
    expect_equal(test_mat[,1], as.numeric(w_mu), tolerance = 0.05)
    expect_equal(test_mat[,3], as.numeric(a$skewness), tolerance = 0.05)
    expect_equal(test_mat[,4], as.numeric(a$kurtosis), tolerance = 0.05)
    # test pfft/qfft
    n <- length(gogarch_sample_index)
    test_mat <- matrix(0, ncol = 2, nrow = n)
    for (i in 1:n) {
        px <- pfft(cf, index = i)
        qx <- qfft(cf, index = i)
        test_mat[i,1] <- px(w_mu[i])
        test_mat[i,2] <- qx(test_mat[i,1])
    }
    expect_equal(test_mat[,2], as.numeric(w_mu), tolerance = 0.05)
})

test_that("GOGARCH Filtering [cond_mean]",{
    new_mod <- tsfilter(mod_full, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(as.numeric(tail(new_mod$spec$target$mu, nrow(filtered))), as.numeric(coredata(filtered)))
    expect_equal(NROW(fitted(new_mod)), NROW(mod_full$mu) + NROW(filtered))
    w <- rep(1/6, 6)
    f <- coredata(fitted(new_mod))
    w_mu <- coredata(new_mod$mu) %*% w
    a <- tsaggregate(new_mod, weights = rep(1/6, 6))
    expect_equal(as.numeric(w_mu), as.numeric(a$mu))
    cf <- tsconvolve(new_mod, weights = w, fft_support = c(-1, 1))
    # test dfft
    n <- length(gogarch_sample_index) + length(gogarch_filter_index)
    test_mat <- matrix(0, ncol = 4, nrow = n)
    for (i in 1:n) {
        dx <- dfft(cf, index = i)
        f <- function(x) x * dx(x)
        w_imu <- integrate(f, -1, 1)$value
        test_mat[i,1] <- w_imu
        test_mat[i,2] <- w_mu[i]
        f3 <- function(x) (x - w_imu)^3 * dx(x)
        f4 <- function(x) (x - w_imu)^4 * dx(x)
        fx3 <- integrate(f3, -1, 1, abs.tol = 1e-8)$value
        fx4 <- integrate(f4, -1, 1, abs.tol = 1e-8)$value
        test_mat[i,3] <- fx3/as.numeric(a$sigma[i])^3
        test_mat[i,4] <- fx4/as.numeric(a$sigma[i])^4
    }
    expect_equal(test_mat[,1], as.numeric(w_mu), tolerance = 0.05)
    expect_equal(test_mat[,3], as.numeric(a$skewness), tolerance = 0.05)
    expect_equal(test_mat[,4], as.numeric(a$kurtosis), tolerance = 0.05)
    # test pfft/qfft
    test_mat <- matrix(0, ncol = 2, nrow = n)
    for (i in 1:n) {
        px <- pfft(cf, index = i)
        qx <- qfft(cf, index = i)
        test_mat[i,1] <- px(w_mu[i])
        test_mat[i,2] <- qx(test_mat[i,1])
    }
    expect_equal(test_mat[,2], as.numeric(w_mu), tolerance = 0.05)
})

test_that("GOGARCH Prediction [cond_mean]",{
    p <- predict(mod_full, h = 6, nsim = 50, cond_mean = pred_mean, seed = 100)
    a <- tsaggregate(p, weights = rep(1/6, 6))
    cf <- tsconvolve(p, weights = rep(1/6, 6))
    mu <- skew <- kurt <- matrix(0, ncol = 6, nrow = 500)
    for (i in 1:6) {
        for (j in 1:500)
        {
            dx <- dfft(cf, index = i, sim = j)
            f <- function(x) x * dx(x)
            w_imu <- integrate(f, -1, 1, rel.tol = 1e-8)$value
            mu[j,i] <- w_imu
            f3 <- function(x) (x - w_imu)^3 * dx(x)
            f4 <- function(x) (x - w_imu)^4 * dx(x)
            fx3 <- integrate(f3, -1, 1, abs.tol = 1e-8)$value
            fx4 <- integrate(f4, -1, 1, abs.tol = 1e-8)$value
            skew[j,i] <- fx3/as.numeric(a$sigma[j,i])^3
            kurt[j,i] <- fx4/as.numeric(a$sigma[j,i])^4
        }
    }
    expect_equal(as.numeric(mu), as.numeric(a$mu), tolerance = 0.05)
})

test_that("GOGARCH Simulation [cond_mean]",{

})

test_that("GOGARCH Estimation DR [cond_mean]",{
})
