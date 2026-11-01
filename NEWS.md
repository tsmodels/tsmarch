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
* Documented the first stage variance regressors in the vignettes, including
the choice between the additive and multiplicative (exponential) forms, and
that first stage mean and variance regressors supplied to the forward-looking
methods are treated as a single deterministic path shared across all
simulated paths, so regressor uncertainty is not propagated into the
predictive distribution.
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
* Fixed the sign of the second stage rows of the partitioned hessian and of
the second stage columns of the score matrix in the DCC and Copula GARCH
models, which were computed from the log-likelihood rather than the negative
log-likelihood used for the first stage. As a result `vcov(type = "H")`
returned negative variances for the second stage parameters and `estfun`
returned second stage scores of the wrong sign.
* Fixed the `QMLE` and `NW` sandwich estimators in `vcov` (and hence `summary`
with `vcov_type = "QMLE"` or `"NW"`), which used `B M B` rather than
`B M B'`. Since the partitioned hessian of the two stage estimator is not
symmetric, this understated the second stage standard errors and returned a
non-symmetric matrix. The default `OP` standard errors used by `summary` are
unaffected.
* The second derivatives in the partitioned hessian are now computed by central
differences with steps relative to the parameter magnitude, improving the
accuracy of the cross derivatives with respect to first stage parameters of
small magnitude such as the GARCH intercept. This requires more likelihood
evaluations, so estimation with `return_hessian = TRUE` is somewhat slower.
* The joint likelihood used in the partitioned hessian (and, for the constant
correlation model, in the scores) of the Copula GARCH model now applies the
same truncation of the probability integral transformed values as the
estimation objective, so that the derivatives are those of the function
actually minimised.
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
