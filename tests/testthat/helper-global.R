suppressMessages(suppressWarnings(library(xts)))
suppressMessages(suppressWarnings(library(tsgarch)))
y <- as.xts(dji30retw)
test_series <- c(1,2,3)
sample_index <- 1:1100
test_index <- 1101:1140

x <- lapply(seq_along(test_series), function(i){
    spec <- garch_modelspec(y[sample_index,i], model = "egarch")
    mod <- estimate(spec, keep_tmb = TRUE)
    return(mod)
})
names(x) <- colnames(y)[test_series]
x <- to_multi_estimate(x)

# CGARCH Model
global_cgarch_spec_constant_p <- cgarch_modelspec(x, dynamics = "constant", transformation = "parametric", copula = "mvn", constant_correlation = "spearman")
global_cgarch_spec_constant_spd <- cgarch_modelspec(x, dynamics = "constant", transformation = "spd", copula = "mvn", constant_correlation = "spearman")
global_cgarch_spec_constant_e <- cgarch_modelspec(x, dynamics = "constant", transformation = "empirical", copula = "mvn", constant_correlation = "spearman")

global_cgarch_constant_estimate_p <- estimate(global_cgarch_spec_constant_p, control = list(trace = 0))
global_cgarch_constant_estimate_spd <- estimate(global_cgarch_spec_constant_spd, control = list(trace = 0))
global_cgarch_constant_estimate_e <- estimate(global_cgarch_spec_constant_e, control = list(trace = 0))

global_cgarch_spec_dcc_normal_p <- cgarch_modelspec(x, dynamics = "dcc", transformation = "parametric", copula = "mvn")
global_cgarch_spec_dcc_normal_spd <- cgarch_modelspec(x, dynamics = "dcc", transformation = "spd", copula = "mvn")
global_cgarch_spec_dcc_normal_e <- cgarch_modelspec(x, dynamics = "dcc", transformation = "empirical", copula = "mvn")

global_cgarch_dcc_estimate_p <- estimate(global_cgarch_spec_dcc_normal_p, control = list(trace = 0), return_hessian = FALSE)
global_cgarch_dcc_estimate_spd <- estimate(global_cgarch_spec_dcc_normal_spd, control = list(trace = 0), return_hessian = FALSE)
global_cgarch_dcc_estimate_e <- estimate(global_cgarch_spec_dcc_normal_e, control = list(trace = 0), return_hessian = FALSE)

# DCC Model
global_dcc_spec_constant <- dcc_modelspec(x, dynamics = "constant", distribution = "mvn")
global_dcc_spec_dynamic <- dcc_modelspec(x, dynamics = "dcc", distribution = "mvn")
global_adcc_spec_dynamic <- dcc_modelspec(x, dynamics = "adcc", distribution = "mvn")

global_dcc_constant_estimate <- estimate(global_dcc_spec_constant, control = list(trace = 0))
global_dcc_dynamic_estimate <- estimate(global_dcc_spec_dynamic, control = list(trace = 0))
global_adcc_dynamic_estimate <- estimate(global_adcc_spec_dynamic, control = list(trace = 0))

# GOGARCH model tested in folder longtests
