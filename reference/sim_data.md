# Simulate acre-formatted SCR data.

Simulate acre-formatted SCR data.

## Usage

``` r
sim_data(
  sim_name,
  n.rand,
  fit = NULL,
  seed = 810,
  suppress_messages = FALSE,
  proportion_missing = 0
)
```

## Arguments

- sim_name:

  a string denoting type of data to be simulated. Supported types can be
  found using the `get_dataset_names()` function.

- n.rand:

  a numeric value denoting the number of data sets to be generated. If
  n.rand \> 1, `output$capt` will be a list of length `n.rand`.

- fit:

  an object generated from the model fitting function "fit.acre()" or
  the bootstrap process "boot.acre()". If `fit` is provided, then all
  parameters will be taken from the fitted object.

- seed:

  a numeric value denoting the random seed, defaults to 810.

- suppress_messages:

  a logical value which indicates whether to suppress helper messages.

- proportion_missing:

  a numeric value used to set a proportion of the covariate data
  generated to NA (currently only bearing, distance and toa supported)

## Value

A data frame containing the simulated data, as well as arguments used.
