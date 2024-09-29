test_that("GOGARCH residuals",{
    stdize <- FALSE
    type <- "standard"
    r <- residuals(global_gogarch_mod, standardize = stdize, type = type)
    expect_equal(as.numeric(r), as.numeric(global_gogarch_mod$spec$target$y))
    stdize <- TRUE
    type <- "whitened"
    r <- coredata(residuals(global_gogarch_mod, standardize = stdize, type = type))
    r <- cor(r)
    r <- r[lower.tri(r)]
    expect_lt(mean(abs(r)), 0.02)
})

test_that("GOGARCH ICA",{
    stdize <- FALSE
    type <- "standard"
    r_1 <- coredata(residuals(global_gogarch_mod, standardize = stdize, type = type))
    r_2 <- coredata(residuals(global_gogarch_mod_dr, standardize = stdize, type = type))
    y <- global_gogarch_mod$spec$target$y
    cor_1 <- diag(cor(r_1, y))
    cor_2 <- diag(cor(r_2, y))
    expect_equal(as.numeric(cor_1), rep(1, 6))
    expect_lt(mean(as.numeric(cor_2)), 1)
})


