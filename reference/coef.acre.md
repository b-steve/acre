# Extract coefficients from the output of acre model

Extract coefficients from the output of acre model

## Usage

``` r
# S3 method for class 'acre'
coef(object, types = NULL, pars = NULL, new.covariates = NULL, ...)
```

## Arguments

- object:

  a fitted model from fit.acre().

- types:

  a character vector, accept any subset from ('all', 'fitted', 'linked',
  'derived'), default is 'linked'. In details: linked - There is a link
  function attached to each parameter, either 'log', 'logit' or
  'identical'. The 'linked' estimation is the estimation before
  back-transformation according to the link functions. For any extended
  parameters, the estimations of the coefficients of any covariates will
  be always under link function of identity. For example, if we have D ~
  x1, the estimations of D_int and D_x1 are 1 and 0.5, if the value of
  'x1' is not provided in the argument of 'new.covariates', then no
  matter what the 'types' is, these estimations will be shown all the
  time. If x1 = 3, then when 'linked', the function will return 1 + 0.5
  \* 3 as the estimation of log(D); when 'fitted', the function will
  return exp(1 + 0.5 \* 3) as the estimation of D. fitted - The back
  transformed estimations from the 'linked' estimations and their link
  function. derived - The estimations of 'esa' for each session.

- pars:

  a character vector containing any parameter names.

- new.covariates:

  a data frame containing the values of covariates of any extended
  parameter.

## Value

a named numeric vector
