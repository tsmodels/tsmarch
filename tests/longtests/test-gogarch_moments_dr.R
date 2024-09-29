test_that("GOGARCH Estimation Coskewness",{
    series <- 6
    stdize <- TRUE
    # distribution not an active argument for estimated object
    index <- 1
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,length(index)))
    index <- 10:15
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,length(index)))

    index <- NULL
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series, length(gogarch_sample_index)))

    index <- 1
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,length(index)))
    index <- 10:15
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,length(index)))
    index <- NULL
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2, length(gogarch_sample_index)))


    stdize <- FALSE
    index <- 1
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,length(index)))
    index <- 10:15
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,length(index)))

    index <- NULL
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series, length(gogarch_sample_index)))

    index <- 1
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,length(index)))
    index <- 10:15
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,length(index)))
    index <- NULL
    k <- tscoskew(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2, length(gogarch_sample_index)))
})

test_that("GOGARCH Estimation Cokurtosis",{
    series <- 6
    stdize <- TRUE
    # distribution not an active argument for estimated object
    index <- 1
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,length(index)))
    index <- 10:15
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,length(index)))

    index <- NULL
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,length(gogarch_sample_index)))

    index <- 1
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,length(index)))
    index <- 10:15
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,length(index)))
    index <- NULL
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3, length(gogarch_sample_index)))


    stdize <- FALSE
    index <- 1
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,length(index)))
    index <- 10:15
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,length(index)))

    index <- NULL
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,length(gogarch_sample_index)))

    index <- 1
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,length(index)))
    index <- 10:15
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,length(index)))
    index <- NULL
    k <- tscokurt(global_gogarch_mod_dr, index = index, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3, length(gogarch_sample_index)))
})

test_that("GOGARCH Prediction Coskewness [h=1]",{
    h <- 1
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
    stdize <- FALSE
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
})

test_that("GOGARCH Prediction Coskewness [h>1]",{
    h <- 3
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
    stdize <- FALSE
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
})


test_that("GOGARCH Simulation Coskewness [h=1]",{
    h <- 1
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100, burn = 100)
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
    stdize <- FALSE
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
})

test_that("GOGARCH Simulation Coskewness [h>1]",{
    h <- 3
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100, burn = 100)
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
    stdize <- FALSE
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,h))
    k <- tscoskew(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h,nsim))
    k <- tscoskew(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^2,h))
})


test_that("GOGARCH Prediction Cokurtosis [h=1]",{
    h <- 1
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
    stdize <- FALSE
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
})

test_that("GOGARCH Prediction Cokurtosis [h>1]",{
    h <- 3
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
    stdize <- FALSE
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
})


test_that("GOGARCH Simulation Cokurtosis [h=1]",{
    h <- 1
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
    stdize <- FALSE
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
})

test_that("GOGARCH Simulation Cokurtosis [h>1]",{
    h <- 3
    nsim <- 10
    series <- 6
    stdize <- TRUE
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
    stdize <- FALSE
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = TRUE)
    expect_equal(dim(k), c(series,series,series,series,h))
    k <- tscokurt(p, distribution = TRUE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h,nsim))
    k <- tscokurt(p, distribution = FALSE, standardize = stdize, folded = FALSE)
    expect_equal(dim(k), c(series,series^3,h))
})


test_that("GOGARCH Covariance/Correlation Estimation",{
    series <- 6
    stdize <- TRUE
    V <- tscov(global_gogarch_mod_dr)
    expect_equal(dim(V), c(series,series,length(gogarch_sample_index)))
    C <- tscor(global_gogarch_mod_dr)
    expect_equal(dim(C), c(series,series,length(gogarch_sample_index)))
    diagC <- C[1,1,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
})

test_that("GOGARCH Prediction Covariance/Correlation",{
    series <- 6
    nsim <- 10
    h <- 1
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(series,series,h, nsim))
    # 1 step ahead no uncertainty
    expect_equal(mean(V[1,1,,]), V[1,1,1,1])
    V <- tscov(p, distribution = FALSE)
    expect_equal(dim(V), c(series,series,h))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(series,series,h,nsim))
    diagC <- C[1,1,,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
    C <- tscor(p, distribution = FALSE)
    expect_equal(dim(C), c(series,series,h))
    diagC <- C[1,1,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)

    series <- 6
    nsim <- 10
    h <- 3
    p <- predict(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(series,series,h, nsim))
    # 1 step ahead no uncertainty
    expect_equal(mean(V[1,1,1,]), V[1,1,1,1])
    V <- tscov(p, distribution = FALSE)
    expect_equal(dim(V), c(series,series,h))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(series,series,h,nsim))
    diagC <- C[1,1,,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
    C <- tscor(p, distribution = FALSE)
    expect_equal(dim(C), c(series,series,h))
    diagC <- C[1,1,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
})


test_that("GOGARCH Simulation Covariance/Correlation",{
    series <- 6
    nsim <- 10
    h <- 1
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(series,series,h, nsim))
    # 1 step ahead no uncertainty
    expect_equal(mean(V[1,1,,]), V[1,1,1,1])
    V <- tscov(p, distribution = FALSE)
    expect_equal(dim(V), c(series,series,h))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(series,series,h,nsim))
    diagC <- C[1,1,,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
    C <- tscor(p, distribution = FALSE)
    expect_equal(dim(C), c(series,series,h))
    diagC <- C[1,1,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)

    series <- 6
    nsim <- 10
    h <- 3
    p <- simulate(global_gogarch_mod_dr, h = h, nsim = nsim, seed = 100)
    V <- tscov(p, distribution = TRUE)
    expect_equal(dim(V), c(series,series,h, nsim))
    # 1 step ahead no uncertainty
    expect_equal(mean(V[1,1,1,]), V[1,1,1,1])
    V <- tscov(p, distribution = FALSE)
    expect_equal(dim(V), c(series,series,h))
    C <- tscor(p, distribution = TRUE)
    expect_equal(dim(C), c(series,series,h,nsim))
    diagC <- C[1,1,,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
    C <- tscor(p, distribution = FALSE)
    expect_equal(dim(C), c(series,series,h))
    diagC <- C[1,1,]
    expect_equal(max(diagC), 1)
    expect_equal(min(diagC), 1)
})
