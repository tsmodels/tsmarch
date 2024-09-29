suppressMessages(suppressWarnings(library(xts)))
suppressMessages(suppressWarnings(library(tsgarch)))
suppressMessages(suppressWarnings(library(future)))
RcppParallel::setThreadOptions(numThreads = 1)
plan(multisession, workers = 5)
options(future.globals.maxSize = 850 * 1024^2)
# GOGARCH Model
gogarch_sample_index <- 1:1600
gogarch_filter_index <- 1608:1690
data("globalindices", package = "tsmarch")
globalindices <- as.xts(globalindices)
print("running long tests")
global_gogarch_spec <- gogarch_modelspec(globalindices[gogarch_sample_index,1:6], distribution = "gh", components = 6, lambda_range = c(-3, 3), shape_range = c(0.25, 20))
global_gogarch_mod <- suppressWarnings(estimate(global_gogarch_spec))
# dimensionality reduced
global_gogarch_spec_dr <- gogarch_modelspec(globalindices[gogarch_sample_index,1:6], distribution = "gh", components = 4, lambda_range = c(-3, 3), shape_range = c(0.25, 20))
global_gogarch_mod_dr <- suppressWarnings(estimate(global_gogarch_spec_dr))


arima_model <- lapply(1:6, function(i){
    arima(as.numeric(globalindices[gogarch_sample_index,i]), order = c(6, 0, 0))
})
res <- do.call(cbind, lapply(1:6, function(i){
    as.numeric(residuals(arima_model[[i]]))
}))
colnames(res) <- colnames(globalindices[,1:6])
res <- xts(as.matrix(res), index(globalindices[gogarch_sample_index]))
mu <- globalindices[gogarch_sample_index,1:6] - res
colnames(mu) <- colnames(globalindices[,1:6])
spec_full <- gogarch_modelspec(globalindices[gogarch_sample_index,1:6], distribution = "gh", components = 6, lambda_range = c(-3, 3), shape_range = c(0.25, 20), cond_mean = coredata(mu))
mod_full <- suppressWarnings(estimate(spec_full))

spec_dr <- gogarch_modelspec(globalindices[gogarch_sample_index,1:6], distribution = "gh", components = 4, lambda_range = c(-3, 3), shape_range = c(0.25, 20), cond_mean = coredata(mu))
mod_dr <- suppressWarnings(estimate(spec_dr))

pred_mean <- do.call(cbind, lapply(1:6, function(i){
    ptmp <- predict(arima_model[[i]], n.ahead = 6)
    as.numeric(ptmp$pred)
}))

filt <- lapply(1:6, function(i){
    cfixed <- coef(arima_model[[i]])
    new_mod <- arima(globalindices[c(gogarch_sample_index, gogarch_filter_index),i], order = c(6,0,0), fixed = cfixed, include.mean = TRUE)
    return(new_mod)
})


filter_index <- index(globalindices[gogarch_filter_index])
res_filtered <- do.call(cbind, lapply(1:6, function(i){
    xts(residuals(filt[[i]]), index(globalindices[c(gogarch_sample_index, gogarch_filter_index)]))
}))
colnames(res_filtered) <- colnames(mu)
res_filtered <- res_filtered[filter_index]


filtered <- do.call(cbind, lapply(1:6, function(i){
    xts(as.numeric(globalindices[c(gogarch_sample_index, gogarch_filter_index),i]) - residuals(filt[[i]]), index(globalindices[c(gogarch_sample_index, gogarch_filter_index)]))
}))
colnames(filtered) <- colnames(mu)
filtered <- filtered[filter_index]

# test_dir("tests/longtests/")
