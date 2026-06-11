# Combining all data for plotting and model fitting

Combines detection data, detector locations, and covariate data into a
single object. Performs spatial interpolation for spatial covariates.
The resulting object can be used to plot detection and covariate data
with
[`plot.acre_data()`](https://b-steve.github.io/acre/reference/plot.acre_data.md),
and to fit models with
[`fit.acre()`](https://b-steve.github.io/acre/reference/fit.acre.md).

## Usage

``` r
read.acre(
  captures,
  traps,
  mask = NULL,
  control.mask = list(),
  loc.cov = NULL,
  time.loc.cov = NULL,
  convert.loc2mask = list(),
  session.cov = NULL,
  trap.cov = NULL,
  dist.cov = NULL,
  cue.rates = NULL,
  survey.length = NULL,
  sound.speed = 330
)
```

## Arguments

- captures:

  A data frame with detection data. Further details are available below.

- traps:

  A matrix or a data frame with two columns, or a list of such matrices
  or data frames for a multi-session model. Each row provides the x- and
  y-coordinates of a detector location, in metres. For multi-session
  models, each component of the list provides detector locations for one
  session.

- mask:

  Optional. A mask object, for example created using
  [`create.mask()`](https://b-steve.github.io/acre/reference/create.mask.md).
  If `NULL`, a mask object will be created automatically based on the
  `traps` and `control.mask` arguments.

- control.mask:

  A list specifying arguments for
  [`create.mask()`](https://b-steve.github.io/acre/reference/create.mask.md),
  other than `traps`, that are used to create a mask object. The
  `buffer` should be provided, and, optionally, the `spacing`.

- loc.cov:

  A data frame, or a list of data frames, containing spatial covariates.
  Data frames must contain columns `x` and `y` for x- and y-coordinates,
  in metres, and additional columns for spatial covariates measured at
  these locations. Missing `NA` values are allowed. Spatial
  interpolation is performed to produce covariate values at the mask
  point locations; see the section on spatial covariates below.

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

- dist.cov:

  A list containing locations of features from which distances are
  calculated, and can be used as spatial covariates. Each component must
  be named after a feature, with a data frame containing columns `x` and
  `y`, recording the the location of a feature.

- cue.rates:

  A numeric vector containing the recorded cue rates of individuals,
  collected separately to the SCR survey. See the section on cue rate
  data, below.

- survey.length:

  A numeric vector or a scalar, containing the length of each session.
  If it is a scalar and there are multiple sessions, the value will be
  assigned to all sessions.

- sound.speed:

  A scalar, the speed of sound in metres per second. This argument is
  only used when time-of-arrival data are included in the `captures`
  data frame. Defaults to 330, the approximate speed of sound in air.

## Captures argument

The captures argument must be a data frame, where each row corresponds
to a single detection of a call by a single detector. The following
columns are required:

- `session`: a session identifier.

- `ID`: a call identifier.

- `trap`: a detector identifier.

Additional columns can specify auxiliary data collected about the calls.
Missing values of `NA` are allowed, so that an auxiliary data type might
only be available for a subet of calls.

- `bearing`: the estimated bearing (in radians) from the detector to the
  detected call.

- `dist`: the estimated distance (in metres) from the detector to the
  detected call.

- `ss`: the received signal strength of the detected call.

- `toa`: the time of arrival of the detected call at the detector,
  measured in seconds since the start of the survey.

- `animal_ID`: an animal identifier, if animals can be recognised from
  their calls.

When bearings, distances, signal strengths, or times of arrival are
provided, they are automatically used by
[`fit.acre()`](https://b-steve.github.io/acre/reference/fit.acre.md)
during model fitting via the method described by Borchers et al (2015),
which can considerably improve precision of density estimates.

If `animal_ID` is provided, then calling animal density and call
production rates are separately estimated from the SCR data via the
method described by Stevenson et al (2021). Note that animals are
assumed to remain stationary throughout the survey.

## Spatial covariates

When spatial covariates are used to model spatial variation in density,
they must be available for every mask point location during
model-fitting. However, the resolution of the spatial data available
from its primary source is likely to differ from the resolution of the
mask. Spatial interpolation is therefore required to generate covariate
values at the mask points.

If required, `read.acre()` will automatically carry out this spatial
interpolation. The user simply needs to provide whatever spatial data
are available from its original source via the arguments `loc.cov`,
`time.loc.cov`, or both. Measurements do not need to be available at
points falling in a regular grid or anything like that; covariates
measured at a higgledy-piggledy collection of coordinates are fine. Our
spatial interpolation method works by selecting the closest measurements
available and taking a weighted average, using a distance-based
weighting function so that closer measurements have greater influence.

Be warned that estimated density surfaces are only as good as the
interpolated covariate values. If the covariate data provided via
`loc.cov` and `time.loc.cov` have low spatial resolution, then the
interpolated covariates may not closely approximate reality, leading to
poor inference. Users should inspect plots of interpolated covariates
using [`plot()`](https://rdrr.io/r/graphics/plot.default.html) to assess
their suitability. See
[`plot.acre_data()`](https://b-steve.github.io/acre/reference/plot.acre_data.md).

If a user does not want to use our spatial interpolation method, they
can take full control of the procedure using these steps:

- Create a mask themselves using
  [`create.mask()`](https://b-steve.github.io/acre/reference/create.mask.md).

- Carry out spatial interpolation themselves, generating covariate
  values at the mask points.

- Provide the spatially interpolated data as the `loc.cov` and
  `time.loc.cov` arguments.

## Cue rate data

If animal identifiers are not included in the `captures` argument, then
call density (calls produced per unit area per unit time) rather than
population density (individuals per unit area) is estimated by
[`fit.acre()`](https://b-steve.github.io/acre/reference/fit.acre.md).

If cue rates are separately collected and provided via the `cue.rates`
argument, then population density can also be estimated by
[`fit.acre()`](https://b-steve.github.io/acre/reference/fit.acre.md).
See Stevenson et al (2015) for further details.

The `cue.rates` argument should be a vector, where each element is the
observed cue rate of an individual animal. The cue rates must be
provided in calls per unit time, where the time units are equivalent to
those used by `survey.length`. For example, if the survey was 30 minutes
long, the cue rates must be provided in cues per minute if
`survey.length = 30`, but in cues per hour if `survey.length = 0.5`.

## References

Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas, L., and Marques,
T. A. (2015) A unifying model for capture-recapture and distance
sampling surveys of wildlife populations. *Journal of the American
Statistical Association*, *110*(509), 195–204.

Stevenson, B. C., Borchers, D. L., Altwegg, R., Swift, R. J., Gillespie,
D. M., and Measey, G. J. (2015) A general framework for animal density
estimation from acoustic detections across a fixed microphone array.
*Methods in Ecology and Evolution*, *6*(1), 38–48.

Stevenson, B. C., van Dam-Bates, P., Young, C. K. Y., and Measey, J.
(2021) A spatial capture-recapture model to estimate call rate and
population density from passive acoustic surveys. \*Methods in Ecology
and Evolution, *12*(3), 432–442.

## Examples

``` r
if (FALSE) { # \dontrun{
## Getting some data.
example <- get("simple_hhn")
## A simple modelthat
simple.hhn.fit <- read.acre(capt = example$capt, traps = example$traps, control.mask = list(buffer = 30))

## A simple model with a hazard-rate detection function.
example_hr <- get("simple_hr")
simple.hr.fit <- read.acre(capt = example_hr$capt, traps = example_hr$traps, control.mask = list(buffer = 30),
                           detfn = "hr")
} # }
```
