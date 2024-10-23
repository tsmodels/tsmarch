test_that("DCC Constant Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(dcc_constant_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(dcc_constant_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("DCC Dynamic Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(dcc_dynamic_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(dcc_dynamic_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("DCC Dynamic Predict Test",{
    p <- predict(dcc_dynamic_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})

test_that("DCC Constant Predict Test",{
    p <- predict(dcc_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})


test_that("DCC Dynamic Simulation Test",{
    p <- simulate(dcc_dynamic_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

test_that("DCC Constant Simulation Test",{
    p <- simulate(dcc_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

###
test_that("CGARCH Constant Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(cgarch_constant_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(cgarch_constant_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("CGARCH Dynamic Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(cgarch_dcc_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(cgarch_dcc_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("CGARCH Dynamic Predict Test",{
    p <- predict(cgarch_dcc_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})

test_that("CGARCH Constant Predict Test",{
    p <- predict(cgarch_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})


test_that("CGARCH Dynamic Simulation Test",{
    p <- simulate(cgarch_dcc_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

test_that("CGARCH Constant Simulation Test",{
    p <- simulate(cgarch_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})


test_that("DCC Constant Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(dcc_constant_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(dcc_constant_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("DCC Dynamic Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(dcc_dynamic_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(dcc_dynamic_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("DCC Dynamic Predict Test",{
    p <- predict(dcc_dynamic_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})

test_that("DCC Constant Predict Test",{
    p <- predict(dcc_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})


test_that("DCC Dynamic Simulation Test",{
    p <- simulate(dcc_dynamic_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

test_that("DCC Constant Simulation Test",{
    p <- simulate(dcc_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

###
test_that("CGARCH Constant Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(cgarch_constant_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(cgarch_constant_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("CGARCH Dynamic Filter Test",{
    n_init <- NROW(res)
    n_new <- n_init + NROW(res_filtered)
    new_mod <- tsfilter(cgarch_dcc_mod, y = res_filtered, cond_mean = coredata(filtered))
    expect_equal(NROW(fitted(new_mod)), n_new)
    expect_equal(coredata(as.numeric(filtered)), as.numeric(tail(coredata(fitted(new_mod)), NROW(filtered))))
    V_fit <- tscov(cgarch_dcc_mod)
    V_filt <- tscov(new_mod)
    expect_equal(as.numeric(V_fit[,,gogarch_sample_index]), as.numeric(V_filt[,,gogarch_sample_index]))
})

test_that("CGARCH Dynamic Predict Test",{
    p <- predict(cgarch_dcc_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})

test_that("CGARCH Constant Predict Test",{
    p <- predict(cgarch_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    predicted_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(predicted_mean, pred_mean, tolerance = 1e-4)
})


test_that("CGARCH Dynamic Simulation Test",{
    p <- simulate(cgarch_dcc_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

test_that("CGARCH Constant Simulation Test",{
    p <- simulate(cgarch_constant_mod, h = NROW(pred_mean), cond_mean = pred_mean, nsim = 1000, seed = 100)
    simulated_mean <- sapply(1:n_series, function(i) rowMeans(p$mu[,i,]))
    expect_equal(simulated_mean, pred_mean, tolerance = 1e-4)
})

