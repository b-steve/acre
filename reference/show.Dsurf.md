# Plotting an estimated density surface

Plots density surface estimated by a model fitted with the function
[fit.acre](https://b-steve.github.io/acre/reference/fit.acre.md)

## Usage

``` r
show.Dsurf(
  fit,
  session = NULL,
  show.cv = FALSE,
  new.data = NULL,
  D.cov = NULL,
  xlim = NULL,
  ylim = NULL,
  x.pixels = 50,
  y.pixels = 50,
  zlim = NULL,
  scale = 1,
  plot.contours = FALSE,
  add = FALSE,
  convert.loc2mask = NULL,
  arg.col = 100,
  trap.plot = NULL,
  ...
)
```

## Arguments

- fit:

  an object generated from the model fitting function "fit.acre()" or
  the bootstrap process "boot.acre()".

- session:

  The session with the detector array and invidual(s) to be plotted.
  Ignored if the `newdata` argument is provided.

- show.cv:

  Logical. If true, the CV of the density estimate is plotted rather
  than the estimate itself. At present, this will only work if `newdata`
  is also provided.

- new.data:

  A data frame including new mask points and covariate values, from
  which to estimate and plot density estimates for. This allows, for
  example, estimates to be provided for new regions not included in the
  mask used to fit the model. Two columns, named `x` and `y`, must be
  included, providing the x- and y-coordinates of the new mask points.
  Additional columns must provide the covariates used to fit the model.

- xlim:

  a numeric vector with two elements as the range of x-axis.

- ylim:

  a numeric vector with two elements as the range of y-axis.

- zlim:

  A numeric vector of length 2, giving the range of the density contours

- scale:

  By default, density is in animals per hectare. The plotted values are
  multiplied by this argument, allowing for user-specified units. For
  example, setting `scale = 100` results in densities plotted as animals
  per square kilometre.

- plot.contours:

  Logical, if `TRUE`, contours are plotted.

- add:

  a logical value indicates whether to add the lines into the existing
  plot

- convert.loc2mask:

  A list to control the spatial interpolation method used to compute
  covariate values for mask locations based on data provided in
  `loc.cov` and `time.loc.cov`. See the section on spatial covariates
  below.

- arg.col:

  A numeric value, indicating the number of levels to stretch the color
  over

- ...:

  otlp
