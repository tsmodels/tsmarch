
## no mean
test_that("DCC Dynamic Filter Test no_mean",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(dcc_dynamic_mod_nomean, y = res_filtered)
    expect_equal(NROW(fitted(new_mod)), n_new)
    V_fit <- tscov(dcc_dynamic_mod_nomean)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("DCC Dynamic Predict Test no_mean",{
    p <- predict(dcc_dynamic_mod_nomean, h = 5, nsim = 5000, seed = 100)
    predicted_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(dcc_dynamic_mod_nomean$mu))
    expect_equal(predicted_mean, expected_mean, tolerance = 1e-4)
})

test_that("DCC Constant Predict Test no_mean",{
    p <- predict(dcc_constant_mod_nomean, h = 5, nsim = 5000, seed = 100)
    predicted_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(dcc_constant_mod_nomean$mu))
    expect_equal(predicted_mean, expected_mean, tolerance = 1.5e-4)
})

test_that("DCC Dynamic Simulation Test no_mean",{
    p <- simulate(dcc_dynamic_mod_nomean, h = 5, nsim = 5000, seed = 100)
    simulated_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(dcc_dynamic_mod_nomean$mu))
    expect_equal(simulated_mean, expected_mean, tolerance = 1e-4)
})

test_that("DCC Constant Simulation Test no_mean",{
    p <- simulate(dcc_constant_mod_nomean, h = 5, nsim = 5000, seed = 100)
    simulated_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(dcc_constant_mod_nomean$mu))
    expect_equal(simulated_mean, expected_mean, tolerance = 2e-4)
})

###
test_that("CGARCH Dynamic Filter Test no_mean",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(cgarch_dcc_mod_nomean, y = res_filtered)
    expect_equal(NROW(fitted(new_mod)), n_new)
    V_fit <- tscov(cgarch_dcc_mod_nomean)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("CGARCH Dynamic Predict Test no_mean",{
    p <- predict(cgarch_dcc_mod_nomean, h = 5, nsim = 5000, seed = 100)
    predicted_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(cgarch_dcc_mod_nomean$mu))
    expect_equal(predicted_mean, expected_mean, tolerance = 1e-4)
})

test_that("CGARCH Constant Predict Test no_mean",{
    p <- predict(cgarch_constant_mod_nomean, h = 5, nsim = 5000, seed = 100)
    predicted_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(cgarch_constant_mod_nomean$mu))
    expect_equal(predicted_mean, expected_mean, tolerance = 1.5e-4)
})

test_that("CGARCH Dynamic Simulation Test no_mean",{
    p <- simulate(cgarch_dcc_mod_nomean, h = 5, nsim = 5000, seed = 100)
    simulated_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(cgarch_dcc_mod_nomean$mu))
    expect_equal(simulated_mean, expected_mean, tolerance = 1e-4)
})

test_that("CGARCH Constant Simulation Test no_mean",{
    p <- simulate(cgarch_constant_mod_nomean, h = 5, nsim = 5000, seed = 100)
    simulated_mean <- colMeans(sapply(1:n_series, function(i) rowMeans(p$mu[,i,])))
    expected_mean <- as.numeric(colMeans(cgarch_constant_mod_nomean$mu))
    expect_equal(simulated_mean, expected_mean, tolerance = 2e-4)
})
