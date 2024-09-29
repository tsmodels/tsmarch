test_that("GOGARCH [DR] Filter Test",{
    local_gogarch_mod_dr <- tsfilter(global_gogarch_mod_dr, y = globalindices[gogarch_filter_index, 1:6])
    K1 <- tscokurt(local_gogarch_mod_dr)
    K2 <- tscokurt(global_gogarch_mod_dr)
    expect_equal(K1[,,,,1], K2[,,,,1])
    expect_equal(K1[,,,,1600], K2[,,,,1600])
})

test_that("GOGARCH Filter Test",{
    local_gogarch_mod <- tsfilter(global_gogarch_mod, y = globalindices[gogarch_filter_index, 1:6])
    K1 <- tscokurt(local_gogarch_mod)
    K2 <- tscokurt(global_gogarch_mod)
    expect_equal(K1[,,,,1], K2[,,,,1])
    expect_equal(K1[,,,,1600], K2[,,,,1600])
})


