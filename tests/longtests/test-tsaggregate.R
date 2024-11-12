test_that("GOGARCH tsaggregate",{
    series <- 6
    w <- rep(1/6, 6)
    p <- predict(global_gogarch_mod, h = 10, nsim = 15000, seed = 100)
    agg_d <- tsaggregate(p, weights = w, distribution = TRUE)
    agg_v <- tsaggregate(p, weights = w, distribution = FALSE)
    expect_equal(as.numeric(sqrt(colMeans(agg_d$sigma^2))), as.numeric(agg_v$sigma), tolerance = 0.01)
})

test_that("GOGARCH [DR] tsaggregate",{
    series <- 6
    w <- rep(1/6, 6)
    p <- predict(global_gogarch_mod_dr, h = 10, nsim = 15000, seed = 100)
    agg_d <- tsaggregate(p, weights = w, distribution = TRUE)
    agg_v <- tsaggregate(p, weights = w, distribution = FALSE)
    expect_equal(as.numeric(sqrt(colMeans(agg_d$sigma^2))), as.numeric(agg_v$sigma), tolerance = 0.02)
})


