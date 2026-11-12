#' @useDynLib ZTRecur, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr group_by summarise across all_of
#' @importFrom stats cov sd var rnorm rpois reshape
#' @importFrom graphics axis legend lines matlines matplot par polygon
#' @importFrom grDevices adjustcolor rainbow
"_PACKAGE"