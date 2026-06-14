# Title

Title

## Usage

``` r
# S3 method for class 'acre'
predict(
  object,
  type = "response",
  newdata = NULL,
  se.fit = TRUE,
  confidence = TRUE,
  level = 0.95,
  realnames = NULL,
  ...
)
```

## Arguments

- object:

  a fitted model from "fit.acre()".

- type:

  a character vector, could be either "response" or "link". The
  "response" correspondence to the type of "fitted" in the function
  "coef.acre()" and the "link" correspondence to the type of "linked" in
  the function "coef.acre()".

- newdata:

  a data frame, the same as "coef.acre()". If there is any extended
  parameters, and the corresponding covariates is not provided here, the
  estimation for such parameters will not be shown.

- se.fit:

  a logical value indicates whether to show standard error.

- confidence:

  a logical value indicates whether to show confidence interval.

- level:

  a numeric value indicates the confident level, default is 0.95.

- realnames:

  a character vector containing any parameter names.

- ...:

  For S3 compatibility.
