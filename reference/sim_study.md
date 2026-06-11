# Runs a simulation study

Runs a simulation study

## Usage

``` r
sim_study(
  sim_name,
  n.rand,
  fit = NULL,
  n.cores = NULL,
  save = T,
  plot = T,
  proportion_missing = 0
)
```

## Arguments

- sim_name:

  a string denoting type of data to be simulated. Supported types can be
  found using the `get_dataset_names()` function.

- n.rand:

  a numeric value denoting the number of data sets to be generated and
  fit. Must be \> 1.

- fit:

  an object generated from the model fitting function "fit.acre()" or
  the bootstrap process "boot.acre()". If `fit` is provided, then all
  parameters will be taken from the fitted object.

- n.cores:

  a numeric value denoting the number of cores to use for simulation.
  Defaults to `floor(parallel::detectCores() * 0.80)`

- save:

  a logical value, indicating whether the sim study data should be
  saved. By default creates a "sim_study" folder to save RDS object to.

- plot:

  a logical value, indicating whether or not to plot the simulated
  estimates.

- proportion_missing:

  a numeric value used to set a proportion of the covariate data to NA.
  Will apply to all covariates present.

## Value

A list of length `n.rand` containing fitted coefficients from the
simulation study.
