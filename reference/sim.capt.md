# Simulating SCR data

Simulates SCR capture histories and associated additional information in
the correct format for use with the function
[fit.acre](https://b-steve.github.io/acre/reference/fit.acre.md). If
`fit` is provided then no other arguments are required. Otherwise, at
least `traps`, `mask`, and `pars` are needed.

## Usage

``` r
sim.capt(
  fit,
  detfn,
  param,
  model = NULL,
  traps,
  control.mask = list(),
  session.cov = NULL,
  trap.cov = NULL,
  loc.cov = NULL,
  dist.cov = NULL,
  time.loc.cov = NULL,
  convert.loc2mask = list(),
  survey.length = NULL,
  ss.opts = NULL,
  cue.rates = NULL,
  n.sessions = NULL,
  n.rand = 1,
  random.location = FALSE,
  sound.speed = 331
)
```

## Arguments

- fit:

  an object generated from the model fitting function "fit.acre()" or
  the bootstrap process "boot.acre()". If `fit` is provided, then all
  other parameters will be taken from the fitted object, so no other
  parameters need to be provided.

- detfn:

  A character string specifying the detection function to be used.
  Either `"hn"` (halfnormal), `"hhn"` (hazard halfnormal), `"hr"`
  (hazard rate), `"th"` (threshold), `"lth"` (log-link threshold), or
  `"ss"` (signal strength). If `"ss"` is used, signal strength
  information must be included in `data`. See the section below on
  parameter names for further details.

- param:

  A named list. Component names are parameter names, and each component
  is the value of the associated parameter. A value for the parameter
  `D`, animal density (or call density, if it an acoustic survey) must
  always be provided, along with values for parameters associated with
  the chosen detection function and additional information type(s).

- model:

  A list with named components. Each component name must match a
  parameter name. The component itself must be a
  [formula](https://rdrr.io/r/stats/formula.html) specifying the
  relationship between covariates and the parameter. See the section on
  model specifications below.

- traps:

  A matrix or a data frame with two columns, or a list of such matrices
  or data frames for a multi-session model. Each row provides the x- and
  y-coordinates of a detector location, in metres. For multi-session
  models, each component of the list provides detector locations for one
  session.

- control.mask:

  A list specifying arguments for
  [`create.mask()`](https://b-steve.github.io/acre/reference/create.mask.md),
  other than `traps`, that are used to create a mask object. The
  `buffer` should be provided, and, optionally, the `spacing`.

- session.cov:

  A data frame containing session covariates. It must contain a column
  `session`, and additional columns for the session-level covariates. If
  spatiotemporal covariates are included using `time.loc.cov`, then a
  column `time` must be included, indicating when the session took
  place.

- trap.cov:

  A data frame containing trap covariates. It must contain a column
  `trap` and additional columns for the trap-level covariates. The
  column `session` must be included for multisession data.

- loc.cov:

  A data frame, or a list of data frames, containing spatial covariates.
  Data frames must contain columns `x` and `y` for x- and y-coordinates,
  in metres, and additional columns for spatial covariates measured at
  these locations. Missing `NA` values are allowed. Spatial
  interpolation is performed to produce covariate values at the mask
  point locations; see the section on spatial covariates below.

- dist.cov:

  A list containing locations of features from which distances are
  calculated, and can be used as spatial covariates. Each component must
  be named after a feature, with a data frame containing columns `x` and
  `y`, recording the the location of a feature.

- time.loc.cov:

  A list, or list of data frames, containing spatial covariates that
  change over time. The data frames must contain columns `x`, `y` and
  `time`, and additional columns for the spatiotemporal covariates. When
  these spatiotemporal covariates are provided, the column `time` must
  also appear in `session.cov`.

- convert.loc2mask:

  A list to control the spatial interpolation method used to compute
  covariate values for mask locations based on data provided in
  `loc.cov` and `time.loc.cov`. See the section on spatial covariates
  below.

- survey.length:

  A numeric vector or a scalar, containing the length of each session.
  If it is a scalar and there are multiple sessions, the value will be
  assigned to all sessions.

- ss.opts:

  A list with information required to fit models that include signal
  strengths as auxiliary detection data. One component must be named
  `cutoff`, a detection threshold. An optional component is `ss.link`,
  which specifies the relationship between distance and the expected
  received signal strength. See the section below on signal strength
  models for further details.

- cue.rates:

  A numeric vector containing the recorded cue rates of individuals,
  collected separately to the SCR survey. See the section on cue rate
  data, below.

- n.sessions:

  An integer for the number of sessions.

- sound.speed:

  A scalar, the speed of sound in metres per second. This argument is
  only used when time-of-arrival data are included in the `captures`
  data frame. Defaults to 330, the approximate speed of sound in air.
