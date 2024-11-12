suppressMessages(suppressWarnings(library(xts)))
suppressMessages(suppressWarnings(library(tsgarch)))
suppressMessages(suppressWarnings(library(future)))
RcppParallel::setThreadOptions(numThreads = 2)
plan(multisession, workers = 6)
options(future.globals.maxSize = 850 * 1024^2)
# GOGARCH Model
gogarch_sample_index <- 1:1600
gogarch_filter_index <- 1608:1690
data("globalindices", package = "tsmarch")
globalindices <- as.xts(globalindices)
print("running long tests")
n_series <- 6
global_gogarch_spec <- gogarch_modelspec(globalindices[gogarch_sample_index,1:n_series], distribution = "gh", components = n_series, lambda_range = c(-3, 3), shape_range = c(0.25, 20))
global_gogarch_mod <- suppressWarnings(estimate(global_gogarch_spec))
# dimensionality reduced
global_gogarch_spec_dr <- gogarch_modelspec(globalindices[gogarch_sample_index,1:n_series], distribution = "gh", components = 4, lambda_range = c(-3, 3), shape_range = c(0.25, 20))
global_gogarch_mod_dr <- suppressWarnings(estimate(global_gogarch_spec_dr))


arima_model <- lapply(1:n_series, function(i){
    arima(as.numeric(globalindices[gogarch_sample_index,i]), order = c(6, 0, 0))
})
res <- do.call(cbind, lapply(1:n_series, function(i){
    as.numeric(residuals(arima_model[[i]]))
}))
colnames(res) <- colnames(globalindices[,1:n_series])
res <- xts(as.matrix(res), index(globalindices[gogarch_sample_index]))
mu <- globalindices[gogarch_sample_index,1:n_series] - res
colnames(mu) <- colnames(globalindices[,1:n_series])

pred_mean <- do.call(cbind, lapply(1:n_series, function(i){
    ptmp <- predict(arima_model[[i]], n.ahead = 5)
    as.numeric(ptmp$pred)
}))

filt <- lapply(1:n_series, function(i){
    cfixed <- coef(arima_model[[i]])
    new_mod <- arima(globalindices[c(gogarch_sample_index, gogarch_filter_index),i], order = c(6,0,0), fixed = cfixed, include.mean = TRUE)
    return(new_mod)
})


filter_index <- index(globalindices[gogarch_filter_index])
res_filtered <- do.call(cbind, lapply(1:n_series, function(i){
    xts(residuals(filt[[i]]), index(globalindices[c(gogarch_sample_index, gogarch_filter_index)]))
}))
colnames(res_filtered) <- colnames(mu)
res_filtered <- res_filtered[filter_index]

filtered <- do.call(cbind, lapply(1:n_series, function(i){
    xts(as.numeric(globalindices[c(gogarch_sample_index, gogarch_filter_index),i]) - residuals(filt[[i]]), index(globalindices[c(gogarch_sample_index, gogarch_filter_index)]))
}))
colnames(filtered) <- colnames(mu)
filtered <- filtered[filter_index]


spec_full <- gogarch_modelspec(res[gogarch_sample_index,1:n_series], distribution = "gh", components = n_series, lambda_range = c(-3, 3), shape_range = c(0.25, 20), cond_mean = coredata(mu))
mod_full <- suppressWarnings(estimate(spec_full))

spec_dr <- gogarch_modelspec(res[gogarch_sample_index,1:n_series], distribution = "gh", components = 4, lambda_range = c(-3, 3), shape_range = c(0.25, 20), cond_mean = coredata(mu))
mod_dr <- suppressWarnings(estimate(spec_dr))


univariate_model <- lapply(1:n_series, function(i){
    spec <- garch_modelspec(res[,i], model = "egarch")
    mod <- estimate(spec, keep_tmb = TRUE)
    return(mod)
})
names(univariate_model) <- colnames(res)
univariate_model <- to_multi_estimate(univariate_model)

# CGARCH Model
cgarch_spec_constant <- cgarch_modelspec(univariate_model, dynamics = "constant", transformation = "parametric", copula = "mvt", constant_correlation = "kendall", cond_mean = mu)
cgarch_constant_mod <- estimate(cgarch_spec_constant, control = list(trace = 0))

cgarch_spec_dcc <- cgarch_modelspec(univariate_model, dynamics = "dcc", transformation = "parametric", copula = "mvt", cond_mean = mu)
cgarch_dcc_mod <- estimate(cgarch_spec_dcc, control = list(trace = 0), return_hessian = FALSE)

# DCC Model
dcc_spec_constant <- dcc_modelspec(univariate_model, dynamics = "constant", distribution = "mvt", cond_mean = mu)
dcc_constant_mod <- estimate(dcc_spec_constant, control = list(trace = 0))

dcc_spec_dynamic <- dcc_modelspec(univariate_model, dynamics = "adcc", distribution = "mvt", cond_mean = mu)
dcc_dynamic_mod <- estimate(dcc_spec_dynamic, control = list(trace = 0), return_hessian = FALSE)



# CGARCH Model
cgarch_spec_constant_nomean <- cgarch_modelspec(univariate_model, dynamics = "constant", transformation = "parametric", copula = "mvt", constant_correlation = "kendall")
cgarch_constant_mod_nomean <- estimate(cgarch_spec_constant_nomean, control = list(trace = 0))

cgarch_spec_dcc_nomean <- cgarch_modelspec(univariate_model, dynamics = "dcc", transformation = "parametric", copula = "mvt")
cgarch_dcc_mod_nomean <- estimate(cgarch_spec_dcc_nomean, control = list(trace = 0), return_hessian = FALSE)

# DCC Model
dcc_spec_constant_nomean <- dcc_modelspec(univariate_model, dynamics = "constant", distribution = "mvt")
dcc_constant_mod_nomean <- estimate(dcc_spec_constant_nomean, control = list(trace = 0))

dcc_spec_dynamic_nomean <- dcc_modelspec(univariate_model, dynamics = "adcc", distribution = "mvt")
dcc_dynamic_mod_nomean <- estimate(dcc_spec_dynamic_nomean, control = list(trace = 0), return_hessian = FALSE)
# test_dir("tests/longtests/")
