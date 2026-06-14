# Extract confidence interval for acre.tmb models

Extract confidence interval for acre.tmb models

## Usage

``` r
# S3 method for class 'acre'
confint(
  object,
  parm = NULL,
  level = 0.95,
  types = NULL,
  new.covariates = NULL,
  ...
)
```

## Arguments

- object:

  a fitted model from "fit.acre()".

- parm:

  a character vector, the same as "coef.acre()".

- level:

  a numeric value indicates the confident level, default is 0.95.

- types:

  a character vector, the same as "coef.acre()".

- new.covariates:

  a data frame, the same as "coef.acre()".

- ...:

  For S3 compatibility.

## Value

a matrix
