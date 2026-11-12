# ZTRecur 0.2.0

* Initial release of `ZTRecur` to CRAN.
* Provides functions for marginal regression and nonparametric estimation of zero-truncated recurrent event data with coarsened censoring times.
* Features C++/Armadillo integration via `Rcpp` and `RcppArmadillo` for efficient computation.
* Key exported functions include:
  * `birth.interval()`: Generates birthdate intervals based on integer ages.
  * `nonpara_esti()`: Nonparametric estimation of marginal mean and rate functions.
  * `plot_nonpar()`: Visualization of estimated marginal mean and rate functions.