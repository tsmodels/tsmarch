# tsmarch 1.0.3

* Added support for joint ARMA(X)-GARCH first stage models in the DCC and Copula
GARCH models, with the conditional mean propagated to the joint distribution
(estimation, filtering, prediction and simulation).
* The `cond_mean` argument is now mutually exclusive with first stage ARMA
dynamics at the model level: if any first stage series uses `arma`, supplying
`cond_mean` raises an error in the specification, filter, simulate and predict
methods.
* Added support for first stage mean equation regressors (`xreg` in
`garch_modelspec`) in the DCC and Copula GARCH models: `tsfilter` and
`predict` gain a `newxreg` argument and `simulate` an `xreg` argument, each
taking a list with one element per series (fully named and possibly partial,
or fully unnamed and complete).
* The `cond_mean` argument is now also mutually exclusive with first stage
mean regressors, for the same reason as with ARMA dynamics.
* Added support for first stage variance equation regressors (`vreg` in
`garch_modelspec`) in the DCC and Copula GARCH models: `tsfilter` and
`predict` gain a `newvreg` argument and `simulate` a `vreg` argument, each
taking a list with one element per series of raw regressor matrices with
`NROW(y)`, `h` and `h + burn` rows respectively. Regressors required by a
first stage model but not supplied raise an error naming the affected
series.
* Fixed `simulate` with `burn > 0` failing for the DCC and Copula GARCH
models: the burn-in period is now discarded once at the end of the joint
recursion rather than in the correlation recursion before the first stage
GARCH recursion could consume it, and the returned correlation, covariance
and standardized residual arrays are correctly of length `h`.
* Minimum required version of tsgarch raised to 1.0.5.
* Fixed the univariate pre-sample trim length in the partitioned hessian/score
calculation when the ARMA order exceeds the GARCH order.
* Fixed the joint hessian/score calculation failing when the first stage models
have unequal numbers of estimated parameters (e.g. when variance targeting is
used for only some series) (#7, reported by @kbuchardt).
* Fixed `pca_cov = "LW"` failing in GOGARCH estimation due to a missing `trace`
argument (#6, reported by @Caliani21).
* Replaced `Rf_error` with `Rcpp::stop` in the C++ code (#5, reported by
@Enchufa2).

# tsmarch 1.0.2

* Fix to DCC model simulation when order > 1.
* Bump requirement for R > 4.1.0 (required by tsgarch)
* Fix to Rf_error prompted from Rcpp team

# tsmarch 1.0.1

* Small fixes to RADICAL algorithm.

* Added the FASTICA algorithm as an option to the GOGARCH model since
the RADICAL algorithm still needs work for high dimensional systems
and may be slow in those cases. This is a custom implementation which
more closely follows the original Matlab code 
(https://research.ics.aalto.fi/ica/fastica/) rather than the other 
alternatives in R.

* Added the constant correlation test of Engle and Sheppard and a flextable
method for pretty printing.
