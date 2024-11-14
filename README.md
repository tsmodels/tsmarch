
# tsmarch <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![R-CMD-check](https://github.com/tsmodels/tsgarch/actions/workflows/rcmdcheck.yaml/badge.svg)](https://github.com/tsmodels/tsmarch/actions/workflows/rcmdcheck.yaml)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--11--14-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsmarch)](https://cran.r-project.org/package=tsmarch)

# tsmarch

The `tsmarch` package represents a re-write and re-think of the models
in [rmgarch](https://CRAN.R-project.org/package=rmgarch). It is written
using simpler S3 methods and classes, has a cleaner code base, extensive
documentation and unit tests, provides speed gains by making use of
parallelization in both R (via the `future` package) and in the C++ code
(via `RcppParallel` package), and works with the new univariate GARCH
package [tsgarch](https://CRAN.R-project.org/package=tsgarch).

## Installation

The package can be installed from CRAN (when published) or the
[tsmodels](github.com/tsmodels/) github repo:

``` r
remotes::install_github("tsmodels/tsmarch", dependencies = TRUE)
```
