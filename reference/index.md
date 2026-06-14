# Package index

## Data preparation

Functions for reading, formatting, and preparing acoustic SCR data.

- [`read.acre()`](https://b-steve.github.io/acre/reference/read.acre.md)
  : Combining all data for plotting and model fitting
- [`create.mask()`](https://b-steve.github.io/acre/reference/create.mask.md)
  : Create mask object
- [`location_density()`](https://b-steve.github.io/acre/reference/location_density.md)
  : Calculate location density

## Model fitting

Fit acoustic SCR models and compute core model quantities.

- [`fit.acre()`](https://b-steve.github.io/acre/reference/fit.acre.md) :
  Fitting acoustic SCR models
- [`esa()`](https://b-steve.github.io/acre/reference/esa.md) : Estimated
  effective sampling areas

## Model summaries

Extract fitted values, parameter estimates, uncertainty, and model
diagnostics.

- [`summary(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/summary.acre.md)
  :

  Summarise `acre` Model Fits

- [`coef(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/coef.acre.md)
  : Extract coefficients from the output of acre model

- [`stdEr()`](https://b-steve.github.io/acre/reference/stdEr.md) : Title

- [`stdEr(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/stdEr.acre.md)
  : Extract standard errors of the estimated parameters from acre models

- [`vcov(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/vcov.acre.md)
  : Extract variance covariance matrix of the estimated parameters from
  acre models

- [`confint(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/confint.acre.md)
  : Extract confidence interval for acre.tmb models

- [`logLik(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/logLik.acre.md)
  : Extract (negative) log-likelihood for a fitted acre object

- [`AIC(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/AIC.acre.md)
  : Calculate AIC for acre Models

- [`predict(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/predict.acre.md)
  : Title

## Bootstrap methods

Bootstrap fitted acre models and summarize bootstrap output.

- [`boot.acre()`](https://b-steve.github.io/acre/reference/boot.acre.md)
  :

  Bootstrapping a fitted `acre` model

- [`coef(`*`<acreboot>`*`)`](https://b-steve.github.io/acre/reference/coef.acreboot.md)
  : Title

- [`stdEr(`*`<acreboot>`*`)`](https://b-steve.github.io/acre/reference/stdEr.acreboot.md)
  : Title

- [`vcov(`*`<acreboot>`*`)`](https://b-steve.github.io/acre/reference/vcov.acreboot.md)
  : Title

- [`confint(`*`<acreboot>`*`)`](https://b-steve.github.io/acre/reference/confint.acreboot.md)
  : Title

- [`predict(`*`<acreboot>`*`)`](https://b-steve.github.io/acre/reference/predict.acreboot.md)
  : Title

## Plotting

Plot acre data, masks, fitted models, and diagnostic outputs.

- [`plot(`*`<acre>`*`)`](https://b-steve.github.io/acre/reference/plot.acre.md)
  : Plotting acre model objects
- [`plot(`*`<acre_data>`*`)`](https://b-steve.github.io/acre/reference/plot.acre_data.md)
  : Plotting acre data

## Simulation

Simulate acoustic SCR data and run simulation studies.

- [`sim.capt()`](https://b-steve.github.io/acre/reference/sim.capt.md) :
  Simulating SCR data
- [`sim_data()`](https://b-steve.github.io/acre/reference/sim_data.md) :
  Simulate acre-formatted SCR data.
- [`sim_study()`](https://b-steve.github.io/acre/reference/sim_study.md)
  : Runs a simulation study
