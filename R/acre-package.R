#' `acre`: Acoustic spatial capture-recapture models
#'
#' The `acre` package provides functions to fit spatial
#' capture-recapture models, with a particular focus on analysing data
#' from acoustic surveys.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib acre
NULL

#' @import Rcpp testthat
#' @importFrom graphics abline arrows axis box contour grid hist image legend lines par plot.new plot.window points polygon text title
#' @importFrom methods is
#' @importFrom stats aggregate coef complete.cases confint cor dgamma dist dnorm median pnorm predict printCoefmat qnorm quantile rgamma rnorm rpois runif setNames sd sigma vcov
#' @importFrom utils data flush.console
NULL
