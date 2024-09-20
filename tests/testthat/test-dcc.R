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
