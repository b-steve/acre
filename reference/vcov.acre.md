# Extract variance covariance matrix of the estimated parameters from acre models

Extract variance covariance matrix of the estimated parameters from acre
models

## Usage

``` r
# S3 method for class 'acre'
vcov(
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

  a fitted model from fit.acre().

- types:

  a character vector, the same as "coef.acre()".

- pars:

  a character vector, the same as "coef.acre()".

- new.covariates:

  a data frame, the same as "coef.acre()".

- show_fixed_par:

  a logical value. To control whether to include the fixed parameters in
  the covariance matrix. It is TRUE by default.

## Value

a list with matrices as its elements if multiple 'types'. a matrix if
only one 'types'.
