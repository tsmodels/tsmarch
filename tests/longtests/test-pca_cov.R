# regression test for issue #6: pca_cov = "LW" failed because the trace
# argument was not passed through .pca() to lw_covariance()
set.seed(10)
n_pca <- 500
y_pca <- xts(matrix(rnorm(n_pca * 3) * 0.01, n_pca, 3), order.by = as.Date("2000-01-01") + seq_len(n_pca))
colnames(y_pca) <- c("A","B","C")

test_that("gogarch estimate works for all pca_cov methods",{
    for (ica in c("radical", "fastica")) {
        for (pc in c("ML", "LW", "EWMA")) {
            spec <- gogarch_modelspec(y_pca, distribution = "nig", model = "garch", order = c(1,1), ica = ica)
            mod <- suppressWarnings(estimate(spec, trace = FALSE, seed = 1000, pca_cov = pc))
            expect_s3_class(mod, "gogarch.estimate")
        }
    }
})
