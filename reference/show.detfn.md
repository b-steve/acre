# Plot the detection function

Plot the detection function

## Usage

``` r
show.detfn(
  fit,
  newdata = NULL,
  skip.extend.param = NULL,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  col = NULL,
  add = FALSE,
  ...
)
```

## Arguments

- fit:

  an object generated from the model fitting function
  [fit.acre](https://b-steve.github.io/acre/reference/fit.acre.md) or
  the bootstrap process
  [boot.acre](https://b-steve.github.io/acre/reference/boot.acre.md).

- newdata:

  data.frame; contains any covariates that will be used for all extended
  parameters (if not be skipped)

- skip.extend.param:

  character; skip extended parameter, for skipped extended parameters,
  use its intercept as the value for this parameter

- xlim:

  a numeric vector with two elements as the range of x-axis.

- ylim:

  a numeric vector with two elements as the range of y-axis.

- main:

  a string as the main title of the plot.

- xlab:

  a string as the sub-title of x-axis.

- ylab:

  a string as the sub-title of y-axis.

- col:

  a string or a numeric vector indicates the color of the plotted line.

- add:

  a logical value indicates whether to add the lines into the existing
  plot.

- ...:

  otlp
