% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{confint.acre_boot}
\alias{confint.acre_boot}
\title{Title}
\usage{
\method{confint}{acre_boot}(
  object,
  types = NULL,
  level = 0.95,
  pars = NULL,
  new.covariates = NULL,
  correct_bias = FALSE,
  from_boot = TRUE,
  ...
)
}
\arguments{
\item{object}{an object generated by the bootstrap function "boot.acre()".}

\item{types}{a character vector, the same as "coef.acre_tmb()".}

\item{level}{a numeric value indicates the confident level, default is 0.95.}

\item{pars}{a character vector, the same as "coef.acre_tmb()".}

\item{new.covariates}{a data frame, the same as "coef.acre_tmb()".}

\item{correct_bias}{a logical value indicates whether to apply a bias correction method to
the confidence interval. Default is FALSE, and the naive percentiles from
the bootstrap results will used as confidence interval; if TRUE, the
bootstrap results will be corrected by beta_boot = 2 * beta_est - beta_boot,
where beta_boot is the bootstrap results values and beta_est is the
estimations from "coef.acre_tmb()", and the percentiles from these corrected
bootstrap results are used as the confidence interval.}

\item{from_boot}{a logical value indicates whether to use the bootstrap results to construct
the confidence matrix. Default is TRUE; if FALSE, "confint.acre_tmb()" will
be called.}

\item{...}{}
}
\description{
Title
}
