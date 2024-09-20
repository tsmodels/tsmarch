
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

test_that("cgarch constant prediction",{

})

test_that("cgarch dynamic prediction",{

})
