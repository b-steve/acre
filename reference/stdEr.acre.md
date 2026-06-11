# Extract standard errors of the estimated parameters from acre models

Extract standard errors of the estimated parameters from acre models

## Usage

``` r
# S3 method for class 'acre'
stdEr(
  object,
  types = NULL,
  pars = NULL,
  new.covariates = NULL,
  show_fixed_par = TRUE,
  ...
)
```

## Arguments

- object:

  a fitted model from "fit.acre()".

- types:

  a character vector, the same as "coef.acre()".

- pars:

  a character vector, the same as "coef.acre()".

- new.covariates:

  a data frame, the same as "coef.acre()".

- show_fixed_par:

  a logical value, the same as "vcov.acre()".

## Value

a named numeric vector
