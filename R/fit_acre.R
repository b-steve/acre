#' Combining all data for plotting and model fitting
#'
#' Combines detection data, detector locations, and covariate data
#' into a single object. Performs spatial interpolation for spatial
#' covariates. The resulting object can be used to plot detection and
#' covariate data with [plot.acre.data()], and to fit models with
#' [fit.acre()].
#'
#' @section Captures argument:
#'
#' The captures argument must be a data frame, where each row
#' corresponds to a single detection, usually of a call, by a single
#' detector. The following columns are required:
#' * `session`: a session identifier.
#' * `ID`: a call identifier.
#' * `trap`: a detector identifier.
#'
#' Additional columns can specify auxiliary data collected about the
#' calls. Missing values of `NA` are allowed, so that an auxiliary
#' data type might only be available for a subet of calls.
#' * `bearing`: the estimated bearing (in radians) from the detector
#'   to the detected call.
#' * `dist`: the estimated distance (in metres) from the detector to
#'  the detected call.
#' * `ss`: the received signal strength of the detected call.
#' * `toa`: the time of arrival of the detected call at the detector,
#'  measured in seconds since the start of the survey.
#' * `animal_ID`: an animal identifier, if animals can be recognised
#' from their calls.
#'
#' When bearings, distances, signal strengths, or times of arrival are
#' provided, they are automatically used by [fit.acre()] during model
#' fitting via the method described by Borchers et al (2015), which
#' can considerably improve precision of density estimates.
#'
#' If `animal_ID` is provided, then calling animal density and call
#' production rates are separately estimated from the SCR data via the
#' method described by Stevenson et al (2021). Note that animals are
#' assumed to remain stationary throughout the survey.
#'
#' @section Spatial covariates:
#'
#'  When spatial covariates are used to model spatial variation in
#'  density, they must be available for every mask point location
#'  during model-fitting. However, the resolution of the spatial data
#'  available from its primary source is likely to differ from the
#'  resolution of the mask. Spatial interpolation is therefore
#'  required to generate covariate values at the mask points.
#'
#'  If required, `read.acre()` will automatically carry out this
#'  spatial interpolation. The user simply needs to provide whatever
#'  spatial data are available from its original source via the
#'  arguments `loc.cov`, `time.loc.cov`, or both. Measurements do not
#'  need to be available at points falling in a regular grid or
#'  anything like that; covariates measured at a higgledy-piggledy
#'  collection of coordinates are fine. Our spatial interpolation
#'  method works by selecting the closest measurements available and
#'  taking a weighted average, using a distance-based weighting
#'  function so that closer measurements have greater influence.
#'
#'  Be warned that estimated density surfaces are only as good as the
#'  interpolated covariate values. If the covariate data provided via
#'  `loc.cov` and `time.loc.cov` have low spatial resolution, then the
#'  interpolated covariates may not closely approximate reality,
#'  leading to poor inference. Users should inspect plots of
#'  interpolated covariates using `plot()` to assess their
#'  suitability. See [plot.acre.data()].
#'
#'  If a user does not want to use our spatial interpolation method,
#'  they can take full control of the procedure using these steps:
#'  * Create a mask themselves using [create.mask()].
#'  * Carry out spatial interpolation themselves, generating covariate
#'    values at the mask points.
#'  * Provide the spatially interpolated data as the `loc.cov` and
#'    `time.loc.cov` arguments.
#'
#' @section Cue rate data:
#'
#' If animal identifiers are not included in the `captures` argument,
#' then call density (calls produced per unit area per unit time)
#' rather than population density (individuals per unit area) is
#' estimated by [fit.acre()].
#'
#' If cue rates are separately collected and provided via the
#' `cue.rates` argument, then population density can also be estimated
#' by [fit.acre()]. See Stevenson et al (2015) for further
#' details.
#'
#' The `cue.rates` argument should be a vector, where each element is
#' the observed cue rate of an individual animal. The cue rates must
#' be provided in calls per unit time, where the time units are
#' equivalent to those used by `survey.length`. For example, if the
#' survey was 30 minutes long, the cue rates must be provided in cues
#' per minute if `survey.length = 30`, but in cues per hour if
#' `survey.length = 0.5`.
#'
#' @references Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas,
#'   L., and Marques, T. A. (2015) A unifying model for
#'   capture-recapture and distance sampling surveys of wildlife
#'   populations. *Journal of the American Statistical Association*,
#'   *110*(509), 195--204.
#' @references Stevenson, B. C., Borchers, D. L., Altwegg, R., Swift,
#'   R. J., Gillespie, D. M., and Measey, G. J. (2015) A general
#'   framework for animal density estimation from acoustic detections
#'   across a fixed microphone array. *Methods in Ecology and Evolution*,
#'   *6*(1), 38--48.
#' @references Stevenson, B. C., van Dam-Bates, P., Young, C. K. Y.,
#'   and Measey, J. (2021) A spatial capture-recapture model to
#'   estimate call rate and population density from passive acoustic
#'   surveys. *Methods in Ecology and Evolution*, *12*(3), 432--442.
#'
#' @param captures A data frame with detection data. Further details
#'   are available below.
#' @param traps A matrix or a data frame with two columns, or a list
#'   of such matrices or data frames for a multi-session model. Each
#'   row provides the x- and y-coordinates of a detector location, in
#'   metres. For multi-session models, each component of the list
#'   provides detector locations for one session.
#' @param mask Optional. A mask object, for example created using
#'   [create.mask()]. If `NULL`, a mask object will be created
#'   automatically based on the `traps` and `control.mask` arguments.
#' @param control.mask A list specifying arguments for
#'   [create.mask()], other than `traps`, that are used to create a
#'   mask object. The `buffer` should be provided, and, optionally,
#'   the `spacing`.
#' @param loc.cov A data frame, or a list of data frames, containing
#'   spatial covariates. Data frames must contain columns `x` and `y`
#'   for x- and y-coordinates, in metres, and additional columns for
#'   spatial covariates measured at these locations. Missing `NA`
#'   values are allowed. Spatial interpolation is performed to produce
#'   covariate values at the mask point locations; see the section on
#'   spatial covariates below.
#' @param time.loc.cov A list, or list of data frames, containing
#'   spatial covariates that change over time. The data frames must
#'   contain columns `x`, `y` and `time`, and additional columns for
#'   the spatiotemporal covariates.  When these spatiotemporal
#'   covariates are provided, the column `time` must also appear in
#'   `session.cov`.
#' @param convert.loc2mask A list to control the spatial interpolation
#'   method used to compute covariate values for mask locations based
#'   on data provided in `loc.cov` and `time.loc.cov`. See the section
#'   on spatial covariates below.
#' @param session.cov A data frame containing session covariates. It
#'   must contain a column `session`, and additional columns for the
#'   session-level covariates.  If spatiotemporal covariates are
#'   included using `time.loc.cov`, then a column `time` must be
#'   included, indicating when the session took place.
#' @param trap.cov A data frame containing trap covariates. It must
#'   contain a column `trap` and additional columns for the trap-level
#'   covariates. The column `session` must be included for
#'   multisession data.
#' @param dist.cov A list containing locations of features from which
#'   distances are calculated, and can be used as spatial
#'   covariates. Each component must be named after a feature, with a
#'   data frame containing columns `x` and `y`, recording the the
#'   location of a feature.
#' @param cue.rates A numeric vector containing the recorded cue rates
#'   of individuals, collected separately to the SCR survey. See the
#'   section on cue rate data, below.
#' @param survey.length A numeric vector or a scalar, containing the
#'   length of each session. If it is a scalar and there are multiple
#'   sessions, the value will be assigned to all sessions.
#' @param sound.speed A scalar, the speed of sound in metres per
#'   second. This argument is only used when time-of-arrival data are
#'   included in the `captures` data frame. Defaults to 330, the
#'   approximate speed of sound in air.
#'
#' @return
#' @export
#'
read.acre <- function(captures, traps, mask = NULL,
                     control.mask = list(), loc.cov = NULL, time.loc.cov = NULL,
                     convert.loc2mask = list(), session.cov = NULL, trap.cov = NULL, dist.cov = NULL,
                     cue.rates = NULL, survey.length = NULL, sound.speed = 330){
  ## Keep all original input arguments.
  arg.names <- names(as.list(environment()))
  arg.input <- vector('list', length(arg.names))
  names(arg.input) <- arg.names
  for(i in arg.names) {
    arg.input[[i]] <- get(i)
  }

  output <- list()
  capt <- create.capt(captures, traps)

  ## The number of sessions: from 'traps' when it has a set of detectors for each session, otherwise
  ## from the captures, as in create.capt().
  if(is(traps, 'list') && !is.data.frame(traps) && length(traps) > 1){
    n.sessions <- length(traps)
  } else {
    n.sessions <- max(capt$session)
  }

  traps <- df.to.list(traps, n.sessions)

  if(!is.null(mask)){
    mask <- df.to.list(mask, n.sessions)

    for(i in 1:n.sessions){
      mask[[i]] <- convert.one.mask(mask[[i]], traps[[i]])
    }

  } else {
    stopifnot(!is.null(control.mask$buffer))
    control.mask$traps <- traps
    mask <- do.call('create.mask', control.mask)
  }



  ## Set up the data for par.extend.
  tem <- par.extend.create(loc.cov = loc.cov, dist.cov = dist.cov, mask = mask,
                                 convert.loc2mask = convert.loc2mask,
                                 session.cov = session.cov, trap.cov = trap.cov,
                                 time.loc.cov = time.loc.cov)

  par.extend <- tem$output
  arg.input$loc.cov <- tem$loc.cov

  output$capt <- capt
  output$traps <- traps
  output$mask <- mask

  output$par.extend <- par.extend
  output$cue.rates <- cue.rates
  output$survey.length <- survey.length
  output$sound.speed <- sound.speed

  ## The object arg.input is not used by fit.from.settings(), but it needs to be passed to the final model fit object.
  output$arg.input <- arg.input
  class(output) <- "acre.data"
  return(output)
}


#' Fitting a spatial capture-recapture model
#'
#' Fits an a spatial capture-recapture model to a data object returned by [read.acre()].
#'
#' @section Model specification:
#'
#' Specify models using the `model` argument.
#'
#' (More to put in here)
#'
#' @section Parameter names:
#'
#' Parameters fall into three groups:
#' * `D`, density (in calls or animals per hectare).
#' * Detection function parameters.
#' * Auxiliary data measurement error parameters.
#'
#' (More to put in here)
#'
#' @section Optimiser settings:
#'
#' The [fit.acre()] function fits a spatial capture-recapture model by
#' numerically maximising the likelihood function with respect to the
#' model parameters. In almost all cases, the optimisation procedure
#' will run without a problem, but in rare cases it might help to
#' adjust how the optimiser works.
#'
#' The `optim.opts` argument allows the user to change some optimiser
#' settings. The argument must be a list. Two optional components
#' control parameter-specific optimisation settings:
#'
#' * `sv`: Overrides default parameter start values.
#' * `fix`: Holds parameters constant rather than estimating them via
#' optimisation.
#'
#' These must be named lists, where component names are parameter
#' names. For `sv`, the component is a scalar that provides the
#' parameter's start value, and for `fix` the component is a scalar
#' that provides the fixed value of the parameter. So, for example,
#' `sv = list(sigma = 200)` will set a start value of 200 m for
#' `sigma`.
#'
#' It is not possible to control coefficients for covariates. These
#' always have start values of zero (i.e., no effect of the
#' covariate), and cannot be fixed. The `sv` component can be named
#' after a parameter that is modelled with covariates. For example, if
#' we include `model = list(D = ~ elevation)`, where `elevation` is
#' height above sea level in metres, then using `sv = list(D = 100)`
#' sets the starting values so that density is spatially homogenous at
#' 100 calls per hectare (or 100 animals per hectare, depending on the
#' model). The `fix` component cannot be used for parameters that are
#' modeled with covariates.
#'
#' Three additional components control optimiser behaviour:
#'
#' * `gr.skip`: A logical value, and `FALSE` by default. We use
#'   automatic differentiation to speed up numeric optimisation, but
#'   in some cases this requires large amounts of RAM. Setting
#'   `gr.skip` to `TRUE` turns off automatic differentiation, which
#'   reduces RAM demands, but, if you have enough RAM, will increase
#'   computation time.
#'
#' * `scale.covs`: A logical value, and `TRUE` by default. If `TRUE`,
#'   covariates are centred and standardised (i.e., we subtract the
#'   mean, then divide by the standard deviation) prior to
#'   optimisation. This leads to an equivalent model, but optimisation
#'   is typically more stable. Regardless of whether `scale.covs` is
#'   `TRUE` or `FALSE`, we report estimated coefficients and standard
#'   errors on the scale of the original, untransformed covariate.
#'
#' * `local`: A logical value, and `FALSE` by default, that toggles
#'   local integration over activity centres. If `TRUE`, the
#'   likelihood will be computed by only integrating over mask points
#'   that are within the buffer distance of all traps that made a
#'   detection. This can speed up optimisation, but makes estimates
#'   more sensitive to the mask buffer setting.
#'
#' @section Signal strength models:
#'
#' There are a few options available for signal strength models.
#'
#' (More to put in here)
#'
#' @param data A data object returned by [read.acre()].
#' @param model A list with named components. Each component name must
#'   match a parameter name. The component itself must be a [formula]
#'   specifying the relationship between covariates and the
#'   parameter. See the section on model specifications below.
#' @param detfn A character string specifying the detection function
#'   to be used. Either `"hn"` (halfnormal), `"hhn"` (hazard
#'   halfnormal), `"hr"` (hazard rate), `"th"` (threshold), `"lth"`
#'   (log-link threshold), or `"ss"` (signal strength).  If `"ss"` is
#'   used, signal strength information must be included in `data`. See
#'   the section below on parameter names for further details.
#' @param ss.opts A list with information required to fit models that
#'   include signal strengths as auxiliary detection data. One
#'   component must be named `cutoff`, a detection threshold. An
#'   optional component is `ss.link`, which specifies the relationship
#'   between distance and the expected received signal strength. See
#'   the section below on signal strength models for further details.
#' @param control.mask a list with elements corresponding to the
#'   function `create.mask`. If provided, and the argument "mask" is
#'   not provided, it will be used together with "traps" in the "dat"
#'   to expand masks.
#' @param mask an object generated by the function `create.mask()`, if
#'   not provided, this function will create it based on provided
#'   `traps` and `control.mask`.
#' @param convert.loc2mask a list. When either the arguments of `mask`
#'   or `control.mask` is provided, the mask will be created by them,
#'   and if there is any parameter extended under mask level, this
#'   argument could be used to control the process of converting
#'   location related covariates to the new mask level data.  For any
#'   details, could refer to the help document in the function
#'   `read.acre()`.
#' @param tracing A logical value. If `TRUE`, information about the
#'   optimiser's journey through the parameter space is printed during
#'   model fitting.
#' @param optim.opts A list with optimisation options. See the section
#'   below on optimiser settings.
#' @param CL A logical value. If `TRUE`, a conditional likelihood
#'   model is fitted, and so density is not estimated.
#' @param two.stage A logical value. If `TRUE`, computations required
#'   to propagate variance under a two-stage density modelling
#'   approach are carried out.
#'
#' @return
#' @export
fit.acre <- function(data, model = NULL, detfn = NULL, ss.opts = NULL,
                    control.mask = NULL, mask = NULL, convert.loc2mask = list(),
                    tracing = TRUE, optim.opts = NULL,
                    two.stage = FALSE, CL = two.stage){
  ## Renaming object.
  dat <- data

  arg.input <- dat$arg.input
  dat$arg.input <- NULL
  mask_override <- FALSE
  ## Pulling out optimiser settings.
  sv <- optim.opts$sv
  bounds <- NULL
  gr.skip <- optim.opts$gr.skip
  fix <- optim.opts$fix
  scale.covs <- optim.opts$scale.covs
  local <- optim.opts$local
  if (is.null(local)){
    local <- FALSE
  }
  if(!is.logical(two.stage) || length(two.stage) != 1L || is.na(two.stage)){
    stop("'two.stage' must be a single TRUE/FALSE value.")
  }

  if(!is.logical(CL) || length(CL) != 1L || is.na(CL)){
    stop("'CL' must be a single TRUE/FALSE value.")
  }

  if(isTRUE(two.stage) && isFALSE(CL)) {
    warning("'two.stage' requires CL = TRUE; setting CL = TRUE.")
    CL <- TRUE
  }

  if(!is.null(mask)){
    dat$mask <- mask
    if(!is.null(control.mask)){
      warning("Argument of 'mask' is provided, 'control.mask' will be ignored.")
    }
    mask_override <- TRUE

  } else if(!is.null(control.mask)){
    stopifnot(!is.null(control.mask$buffer))
    control.mask$traps <- dat$traps
    mask <- do.call("create.mask", control.mask)
    dat$mask <- mask
    mask_override <- TRUE
  }

  if(!is.null(model)){
    if(is.null(dat$par.extend)){
      dat$par.extend <- list(model = model)
    } else {
      dat$par.extend$model <- model
    }
    dat$par.extend$scale <- scale.covs
    dat$par.extend$link <- NULL

    if(mask_override && !is.null(dat$par.extend$data$mask)){
      convert.loc2mask$loc.cov <- arg.input$loc.cov
      convert.loc2mask$mask <- mask
      dat$par.extend$data$mask <- do.call("location.cov.to.mask", convert.loc2mask)
    }

  } else {
    dat$par.extend <- NULL
  }

                                        ## Conditional likelihood checks.
  if(CL){

    if(!is.null(model) && "D" %in% names(model)){
      stop("When CL = TRUE, a density model for 'D' cannot be specified.")
    }

    if(!is.null(dat$par.extend$model) && "D" %in% names(dat$par.extend$model)){
      stop("When CL = TRUE, 'D' cannot appear in the model specification.")
    }

    if(is.null(fix)) fix <- list()

    if(is.null(fix$D)){
      if(!is.null(sv) && !is.null(sv$D)){
        fix$D <- sv$D
      } else {
        fix$D <- 1
      }
    }

    if(!is.null(bounds) && "D" %in% names(bounds)){
      bounds$D <- NULL
    }
  }

  dat$local <- local
  dat$tracing <- tracing
  dat$gr.skip <- gr.skip
  dat$detfn <- detfn
  dat$sv <- sv
  dat$bounds <- bounds
  dat$fix <- fix
  dat$ss.opts <- ss.opts
  sv.link <- NULL
  dat$sv.link <- sv.link
  dat$CL <- CL
  dat$two.stage <- two.stage

  output <- do.call("fit.from.settings", dat)
  output$arg_input <- arg.input
  output$call <- match.call()
  output$CL <- CL

  return(output)
}

#' Fit a model from a list of settings
#'
#' The model-fitting engine used by `fit.acre()`, `boot.acre()` and `sim.study()`. It takes the
#' settings that `fit.acre()` builds from its arguments, which are also stored in a fitted object
#' as `fit$args`, so a model can be refitted to new capture data. Fitting has three stages:
#' `prepare.fit()`, `run.tmb()` and `build.fit()`.
#'
#' @param capt a list or a data frame with capture data. Columns named `session`, `ID` and `trap` are required, with
#'                 each row being regarded as one detection. Extra information can be provided as columns `bearing`, `dist`
#'                 ,`ss`, `toa`.  `animal_ID` could be provided if individuals could be distinguished from their acoustic detection.
#'                 In addition, any individuals from different sessions will be treated as independent ones.
#'                 The default unit for any distance related value should be 'meter', and for time related value should be 'seconds'..
#' @param traps a data frame or a list, the same as the argument 'trap' in function `create.capt()`, please refer to the description
#'              for that function.
#' @param mask a list or a data frame. Generated by the function `create.mask()`.
#' @param detfn same as the argument with the same name in the function `read.acre()`.
#' @param sv same as the argument with the same name in the function `read.acre()`.
#' @param bounds same as the argument with the same name in the function `read.acre()`.
#' @param fix same as the argument with the same name in the function `read.acre()`.
#' @param ss.opts same as the argument with the same name in the function `read.acre()`.
#' @param cue.rates same as the argument with the same name in the function `read.acre()`.
#' @param survey.length same as the argument with the same name in the function `read.acre()`.
#' @param sound.speed same as the argument with the same name in the function `read.acre()`.
#' @param local same as the argument with the same name in the function `fit.acre()`.
#' @param par.extend a list. contains elements of `model`, `data`, `scale` and `link`.
#'                   model: the same as argument `model`  in the function `fit.acre()`;
#'                   data: a list, could contain elements of `session`, `trap`, `mask`. Each of
#'                         them should be a data frame contains the corresponding information
#'                         of sessions, traps or masks and the corresponding covariates;
#'                   scale: the same as argument `is.scale`  in the function `fit.acre()`;
#'                   link: the same as argument `model.link`  in the function `fit.acre()`.
#' @param tracing same as the argument with the same name in the function `fit.acre()`.
#' @param gr.skip same as the argument with the same name in the function `fit.acre()`.
#' @param sv.link a list; this is mostly for development purpose, not recommended to use.
#' @param CL same as the argument with the same name in the function `fit.acre()`.
#' @param two.stage same as the argument with the same name in the function `fit.acre()`.
#'
#' @param ...
#' @keywords internal
fit.from.settings <- function(capt, traps, mask, detfn = NULL, sv = NULL, bounds = NULL, fix = NULL, ss.opts = NULL,
                              cue.rates = NULL, survey.length = NULL, sound.speed = 331, local = FALSE, par.extend = NULL,
                              tracing = TRUE, gr.skip = FALSE, sv.link = NULL, CL = FALSE, two.stage = FALSE, ...){
  ## Keep all original input arguments.
  arg.names <- names(as.list(environment()))
  arg.input <- vector('list', length(arg.names))
  names(arg.input) <- arg.names
  for(i in arg.names) {
    if(!is.null(get(i))){
      arg.input[[i]] <- get(i)
    }
  }

  spec <- prepare.fit(capt = capt, traps = traps, mask = mask, detfn = detfn, sv = sv, bounds = bounds, fix = fix,
                      ss.opts = ss.opts, cue.rates = cue.rates, survey.length = survey.length,
                      sound.speed = sound.speed, local = local, par.extend = par.extend, sv.link = sv.link, CL = CL,
                      arg.input = arg.input, extra.args = list(...))
  tmb <- run.tmb(spec, tracing = tracing, gr.skip = gr.skip, two.stage = two.stage)
  build.fit(spec, tmb)
}

## Prepares a fit: the data, design matrices, start values, bounds and fixed parameters, and the
## data, parameters and map for TMB. Returns everything run.tmb() and build.fit() need.
prepare.fit <- function(capt, traps, mask, detfn, sv, bounds, fix, ss.opts, cue.rates, survey.length, sound.speed,
                        local, par.extend, sv.link, CL, arg.input, extra.args){

  ####################################################################################################
  ## The captures are one row per detection (see create.capt()); with an "animal_ID" column, the
  ## animal model is fitted.
  if(!is.data.frame(capt)) stop('invalid input "capt", "create.capt()" is recommended for obtaining valid input object.')
  animal.model <- 'animal_ID' %in% colnames(capt)
  if(animal.model && !is.null(cue.rates)){
    cue.rates <- NULL
    warning('argument "cue.rates" is ignored since "animal.model" is TRUE.')
  }

  ## The detectors in each session; a single set of detectors is used in every session.
  if(!is(traps, 'list') || is.data.frame(traps)) traps <- list(traps)
  if(length(traps) == 1 && nrow(capt) > 0 && max(capt$session) > 1) traps <- rep(traps, max(capt$session))

  ####################################################################################################
  ## Beginning of fit.acre().
  o <- capture.fun(capt = capt, traps = traps, animal.model = animal.model)
  dims <- o$dims
  data.full <- o$data.capt
  bucket_info <- o$bucket_info

  #########################################################################################################
  data.traps <- trap.fun(traps = traps, dims = dims)
  data.full <- merge(data.full, data.traps, by = c('session', 'trap'), all = TRUE)
  data.full <- data.sort(data.full, 'data.full')

  ########################################################################################################
  o.mask <- mask.fun(mask = mask, dims = dims, animal.model = animal.model,
                    data.traps = data.traps, data.full = data.full,
                    local = local, bucket.info = bucket_info, sound.speed = sound.speed)
  dims[["n.masks"]] <- o.mask$n.masks
  data.full <- o.mask$data.full
  bucket_info <- o.mask$bucket_info
  ## Keep the old version of the output temporarily; it might be useful later.
  mask <- o.mask$mask
  arg.input[['mask']] <- mask

  ## dists = o.mask$dists
  A <- o.mask$A
  buffer <- o.mask$buffer
  ## all.which.local = o.mask$all.which.local
  ## all.n.local = o.mask$all.n.local
  data.dists.thetas <- o.mask$data.dists.thetas
  data.ID_mask <- o.mask$data.ID_mask
  data.mask <- o.mask$data.mask

  ########################################################################################################
  ## In the animal_ID model, the cue rate is no longer an input argument but a fitted parameter, named "mu".
  fulllist.par <- fulllist.par.generator()

  o.par.extend <- par.extend.fun(par.extend = par.extend, data.full = data.full, data.traps = data.traps, data.mask = data.mask,
                                animal.model = animal.model, dims = dims, fulllist.par = fulllist.par,
                                extra.args = extra.args)

  data.full <- o.par.extend$data.full
  data.trap <- o.par.extend$data.trap
  data.mask <- o.par.extend$data.mask
  data.par <- o.par.extend$data.par
  name.extend.par <- o.par.extend$name.extend.par
  gam_output <- o.par.extend$gam_output
  lst_mean_std <- o.par.extend$lst_mean_std
  is.scale <- o.par.extend$is.scale

  ## We didn't convert 'animal_ID' and 'ID' to successive natural numbers earlier because, with 'animal_ID'-level
  ## covariates, the user may want to assign covariates to each 'animal_ID', so we must keep
  ## the original input; but after dealing with 'par.extend', there is no interaction with users about
  ## 'animal_ID' or 'ID', so we convert them here to make sure the design matrices are in the right order.

  ## (Additional note, 09/09/2022: we may not include animal_ID-level covariates, for quite a few reasons;
  ## for example, if we add them, the simulation will not work, as animal_ID is simulated, so the output
  ## animal_ID will be random and we cannot assign any covariate to it. But since the model works,
  ## there is no need to change this part.)

  ## In convert.natural.number(), the key component is as.numeric(as.factor(xxx)). Both 'data.full' and
  ## 'data.ID_mask' contain the same combinations of 'animal_ID' and 'ID', but possibly in a different order.
  ## To make sure as.numeric(as.factor(xxx)) gives the same numeric value for the same 'animal_ID' or 'ID',
  ## sort these two data sets first, so they have the same content in the same order; then there is no
  ## possibility of a different result.
  data.full <- data.sort(data.full, 'data.full')
  data.full <- convert.natural.number(data.full, animal.model, 'both')
  ## Because convert.natural.number() uses rbind(), sort again.
  data.full <- data.sort(data.full, 'data.full')
  if(!is.null(data.ID_mask)){
    data.ID_mask <- data.sort(data.ID_mask, 'data.ID_mask')
    data.ID_mask <- convert.natural.number(data.ID_mask, animal.model, 'both')
    data.ID_mask <- data.sort(data.ID_mask, 'data.ID_mask')
  }

  if(animal.model){
    ## The vector n.animal.call has length sum(n.animals), and each element is the number of calls
    ## made by that animal; since the order of animal_ID matters, we calculate it here.
    tem <- subset(data.full, !is.na(data.full$animal_ID))
    tem <- agg.sort(tem, 'ID', c('session', 'animal_ID'), function(x) length(unique(x)))
    dims$n.animal.call <- tem$x
  } else {
    dims$n.animal.call <- 0
  }


  #####################################################################################################
  o.ss <- ss.fun(ss.opts = ss.opts, data.full = data.full, data.ID.mask = data.ID_mask,
                animal.model = animal.model, dims = dims, bucket.info = bucket_info,
                sv = sv, fix = fix)

  data.full <- o.ss$data.full
  data.ID_mask <- o.ss$data.ID_mask
  dims <- o.ss$dims
  bucket_info <- o.ss$bucket_info
  ss.opts <- o.ss$ss.opts
  ss.link <- ss.opts$ss.link


  ######################################################################################################
  o.CR_SL <- CR.SL(cue.rates = cue.rates, survey.length = survey.length,
                  bucket.info = bucket_info, dims = dims)

  survey.length <- o.CR_SL$survey.length
  ## Update the survey length in the input, so it is included in the output as well.
  arg.input[['survey.length']] <- survey.length

  bucket_info <- o.CR_SL$bucket_info
  ## This is the mean of cue.rate, and it is always generated, but when animal.model = TRUE
  ## this "mean.cue.rates" is not used in the model.
  mean.cue.rates <- o.CR_SL$mean.cue.rates


  ########################################################################################################
  o.param <- param.detfn.fun(animal.model= animal.model, sv = sv, fix = fix, bounds = bounds,
                            name.extend.par = name.extend.par, detfn = detfn, data.full = data.full, data.trap = data.trap,
                            data.mask = data.mask, data.par = data.par, ss.opts = ss.opts,
                            bucket.info = bucket_info, fulllist.par = fulllist.par, A = A, buffer = buffer,
                            survey.length = survey.length, dims = dims, sv.link = sv.link)

  DX.full <- o.param$design.matrices.full
  DX.mask <- o.param$design.matrices.mask
  data.par <- o.param$data.par
  sv.input <- o.param$sv.input
  fix.input <- o.param$fix.input
  bounds.input <- o.param$bounds.input
  name.extend.par <- o.param$name.extend.par
  detfn <- o.param$detfn
  param.og <- o.param$param.og
  ## Make sure param.og follows the order of fulllist.par;
  ## this is important because the "parameter" argument of the TMB model
  ## follows the order of fulllist.par.
  param.og <- fulllist.par[which(fulllist.par %in% param.og)]

  if(is.null(ss.opts)){
    cutoff <- 0
  } else {
    cutoff <- ss.opts$cutoff
  }

  #####################################################################

  ## The dist data are a little special because they could appear in log(capt_dist) in the likelihood.
  ## capt_dist = ifelse(is.na(data.full$dist), 0, data.full$dist)
  ## capt_dist = ifelse(capt_dist == 0, 1e-20, capt_dist)

  if(is.null(data.ID_mask)){
    ## Not used without times of arrival or local integration, but TMB needs the inputs.
    data.ID_mask <- data.frame(local = 1, toa_ssq = 0)
  } else {
    if(!'local' %in% bucket_info){
      data.ID_mask <- data.ID_mask[,-which(colnames(data.ID_mask) == 'local'), drop = FALSE]
    }
    if(!"toa" %in% bucket_info){
      data.ID_mask <- data.ID_mask[,-which(colnames(data.ID_mask) == 'toa_ssq'), drop = FALSE]
    }
  }

  ## Capture data for TMB, one row per call and detector; the placeholder rows that data.full
  ## keeps for sessions without detections are not needed.
  data.calls <- data.full[!is.na(data.full$ID), , drop = FALSE]

  data <- list(n_sessions = dims$n.sessions,
               n_animals = dims$n.animals,
               n_IDs = dims$n.IDs,
               n_traps = dims$n.traps,
               n_masks = dims$n.masks,
               n_calls_each_animal = dims$n.animal.call,


               nrow_data_mask = nrow(data.mask),
               nrow_dx = nrow(data.dists.thetas),

               A = A,
               survey_length = survey.length,
               cue_rates = mean.cue.rates,


               cutoff = cutoff,
               ## Codes for detfn_index:
               ## 1:hn, 2:hhn, 3:hr, 4:th, 5:lth, 6:ss.
               detfn_index = numeric.detfn(detfn),


               ## Simply use the row index of 'data.par' as the numeric substitute for parameters:
               ## 1:g0, 2:sigma, 3:lambda0, 4:z, 5:shape.1, 6:shape.2, 7:shape, 8:scale, 9:b0.ss,
               ## 10:b1.ss, 11:b2.ss, 12:sigma.ss, 13:kappa, 14:alpha, 15:sigma.toa,
               ## 16:sigma.b0.ss, 17:D, 18:mu.
               ## Remember that in TMB indices begin at 0, so all of the indices above should have 1 subtracted.

               param_og = (1:nrow(data.par))[which(data.par[['par']] %in% param.og)],

               par_n_col = as.matrix(data.par[, c('n_col_full', 'n_col_mask')]),
               ## Link index: 1 - identity, 2 - log, 3 - logit.
               par_link = numeric.link(data.par[['link']]),
               ## Link index: 1 - identity, 2 - log, 4 - spherical.
               ss_link = numeric.link(ss.link),

               is_animalID = as.numeric(animal.model),
               is_ss = as.numeric("ss" %in% bucket_info),
               is_bearing = as.numeric("bearing" %in% bucket_info),
               is_toa = as.numeric("toa" %in% bucket_info),
               is_dist = as.numeric("dist" %in% bucket_info),
               is_local = as.numeric("local" %in% bucket_info),


               capt_bin = ifelse(is.na(data.calls$bincapt), NA_real_, data.calls$bincapt),
               capt_bearing = ifelse(is.na(data.calls$bearing), NA_real_, data.calls$bearing),
               capt_dist = ifelse(is.na(data.calls$dist), NA_real_, data.calls$dist),
               capt_ss = ifelse(is.na(data.calls$ss), NA_real_, data.calls$ss),


               dx = data.dists.thetas$dx,
               theta = data.dists.thetas$theta,

               index_local = as.numeric(data.ID_mask$local),
               toa_ssq = as.numeric(data.ID_mask$toa_ssq),

               ## Conditional likelihood.
               is_conditional = as.numeric(CL),
               ## Two-stage fitting (calculate logESA partial derivatives).
               is_two_stage = 0
  )

  ## To avoid "." in the .cpp file.
  fulllist.par.4cpp <- gsub("\\.", "_", fulllist.par)
  param.og.4cpp <- gsub("\\.", "_", param.og)

  parameters <- vector('list', length(fulllist.par) + 1)
  names(parameters) <- c(fulllist.par.4cpp, "u")

  for(i in 1:length(fulllist.par)){
    name.r <- fulllist.par[i]
    name.cpp <- fulllist.par.4cpp[i]
    data[[paste0(name.cpp, "_DX")]] <- ifelse(is.na(DX.full[[name.r]]), 0, DX.full[[name.r]])
    ## The column in the design matrix for the intercept should always be 1.
    data[[paste0(name.cpp, "_DX")]][,1] <- ifelse(data[[paste0(name.cpp, "_DX")]][,1] == 0, 1,
                                                 data[[paste0(name.cpp, "_DX")]][,1])
    data[[paste0(name.cpp, "_DX_mask")]] <- ifelse(is.na(DX.mask[[name.r]]), 0, DX.mask[[name.r]])
    data[[paste0(name.cpp, "_bound")]] <- bounds.input[[name.r]]
    parameters[[name.cpp]] <- sv.input[[name.r]]
  }

  parameters$u <- numeric(sum(dims$n.IDs))

  ## Set the "map" argument for fixed parameters;
  ## obtain the fixed parameters assigned by the user.
  name.fixed.par <- names(fix.input)[!sapply(fix.input, is.null)]
  name.fixed.par.4cpp <- gsub("\\.", "_", name.fixed.par)
  ## TMB cannot optimise with no free parameters.
  if(all(param.og %in% name.fixed.par)){
    stop("All parameters are fixed, so there is nothing to estimate.")
  }

  ## Combine this with the parameters that will not be used in this model.
  name.fixed.par.4cpp <- c(name.fixed.par.4cpp, fulllist.par.4cpp[!(fulllist.par.4cpp %in% param.og.4cpp)])

  ## In case there is any duplication, remove it.
  name.fixed.par.4cpp <- unique(name.fixed.par.4cpp)


  name.fixed.par.4cpp <- c(name.fixed.par.4cpp, "u")


  map <- vector('list', length(name.fixed.par.4cpp))
  names(map) <- name.fixed.par.4cpp

  for(i in 1:length(name.fixed.par.4cpp)){
    par_name <- name.fixed.par.4cpp[i]
    map[[par_name]] <- factor(rep(NA, length(parameters[[par_name]])))
  }

  ## If scaling and centering has been applied to the extended parameters,
  ## we are going to need the back-transformation matrices.
  ## Note that in the case that scaling has not been applied, or a parameter has
  ## not been extended, the associated back transformation matrix is the identity
  ## matrix.
  ## Also note that we need to make sure the parameter names we use match the
  ## output style of TMB.
  TMB_par_names <- rep(gsub("\\.", "_", param.og), lengths(parameters[gsub("\\.", "_", param.og)]))
  G_lst <- construct.G.matrix.list(TMB_par_names, name.extend.par, lst_mean_std,
                                   is.scale)

  list(data = data,
       parameters = parameters,
       map = map,
       G_lst = G_lst,
       name.extend.par = name.extend.par,
       data.trap = data.trap,
       data.mask = data.mask,
       param.og = param.og,
       param.og.4cpp = param.og.4cpp,
       name.fixed.par = name.fixed.par,
       name.fixed.par.4cpp = name.fixed.par.4cpp,
       dims = dims,
       data.par = data.par,
       data.full = data.full,
       data.traps = data.traps,
       data.dists.thetas = data.dists.thetas,
       detfn = detfn,
       DX.full = DX.full,
       DX.mask = DX.mask,
       fix.input = fix.input,
       bucket_info = bucket_info,
       cue.rates = cue.rates,
       mean.cue.rates = mean.cue.rates,
       A = A,
       survey.length = survey.length,
       sound.speed = sound.speed,
       par.extend = par.extend,
       arg.input = arg.input,
       gam_output = gam_output,
       lst_mean_std = lst_mean_std,
       is.scale = is.scale,
       ss.link = ss.link,
       cutoff = cutoff,
       CL = CL)
}

## Builds the TMB objective, maximises the likelihood, and computes standard errors; the estimates
## and their covariance matrix are back-transformed to the original scale of any scaled covariates.
run.tmb <- function(spec, tracing, gr.skip, two.stage){
  data <- spec$data
  parameters <- spec$parameters
  map <- spec$map
  G_lst <- spec$G_lst
  name.extend.par <- spec$name.extend.par
  data.trap <- spec$data.trap
  data.mask <- spec$data.mask
  param.og.4cpp <- spec$param.og.4cpp
  name.fixed.par.4cpp <- spec$name.fixed.par.4cpp
  dims <- spec$dims

  if(!gr.skip) {
    ## Have TMB build the likelihood function.
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
                          silent = TRUE, DLL="acre")
    ## True so that the Hessian is calculated alongside when we optimise.
    obj$hessian <- TRUE

    ## If tracing is enabled, use custom trace output.
    if (tracing) {
      ## If the model is extended, make sure to build the extended parameter names.
      ## These need to correspond to the GLM design matrix column names.
      par_names <- get.coef.names(names(obj$par), name.extend.par,
                                  data.trap, data.mask)

      tap <- make.buffer.printer(trace.cols = par_names,
                                 step = "none",
                                 show.mgc = TRUE)

      ## Capture the originals and wrap them safely.
      orig_fn <- obj$fn
      orig_gr <- obj$gr

      ## Wrap the objective / gradient functions, so each time they are
      ## evaluated, we print the custom trace information.
      obj$fn <- tap.fn(orig_fn, tap, par_names, G_lst)
      obj$gr <- tap.gr(orig_gr, tap, par_names, G_lst)
    }

    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control=list(trace=0))

    ## Restore the originals once the initial optimisation is done
    ## (avoid additional prints when using sdreport()).
    if (tracing) {
      obj$fn <- orig_fn
      obj$gr <- orig_gr
    }

    o <- TMB::sdreport(obj)
  } else {
    ## In the case we don't use TMB to calculate the gradients, we will need to
    ## proceed slightly differently.

    ## The argument type = 'Fun' specifies skipping the gradient calculation.
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
                          silent = TRUE, DLL="acre", type = 'Fun')

    ## Get the names of the fitted parameters,
    ## that is, all the parameters we pass to TMB that are not fixed.
    par_name_fitted <- setdiff(param.og.4cpp, name.fixed.par.4cpp)

    ## The name of ini_par_val is not right; fix it later.
    ini_par_val <- list.2vector.4value(parameters[par_name_fitted])
    fn.base <- function(par) environment(obj$fn)$f(par, type = "double")

    ## If tracing is enabled, use custom trace output.
    if (tracing) {
      ## As before, if the model is extended, make sure to build extended names.
      par_names <- get.coef.names(par_name_fitted, name.extend.par,
                                  data.trap, data.mask)

      tap <- make.buffer.printer(trace.cols = par_names,
                                 step = "max",
                                 show.mgc = FALSE)

      fn_traced <- tap.fn(fn.base, tap, par_names, G_lst)

      opt <- stats::nlminb(ini_par_val, fn_traced, control = list(trace = 0))
    } else {
      opt <- stats::nlminb(ini_par_val, fn.base)
    }

    ## Numerically solve for the Hessian.
    ## Note that the standard errors for g0 may be NaN due to numerical issues,
    ## alongside the fact that its true value lies on the boundary of the
    ## parameter space.
    H <- stats::optimHess(opt$par, fn.base)
    H <- solve(H)
    o <- gr.free.o.restore(fn.base, opt, H, parameters, param.og.4cpp, dims$n.sessions)
  }

  ## If partial derivatives are being returned (i.e. two.stage == T).
  esa_partial_derivatives <- NULL
  if (isTRUE(two.stage)) {
    data$is_two_stage <- as.numeric(TRUE)
    cat("Calculating partial derivatives...\n")
    pd_object <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
                                silent = TRUE, DLL="acre", ADreport = TRUE)
    partial_derivative_matrix <- pd_object$gr(opt$par)
    ## Extracting the names using pd_object$fn() means that in the future, new
    ## ADREPORT()'s won't mess things up (i.e. if we extracted positionally instead).
    log_esa_pd_rows <- which(names(pd_object$fn()) == "log_esa")
    esa_partial_derivatives <- partial_derivative_matrix[log_esa_pd_rows,]
  }

  tmb_output_og <- list(est = o$value, vcov = o$cov)

  ## Restore the estimated parameters' estimates and variance matrix to the original scale.
  G_lst[["esa"]] <- diag(dims$n.sessions)
  est_names <- names(o$value)
  u_names <- unique(est_names)

  missing_G <- setdiff(u_names, names(G_lst))
  for(i in missing_G){
    G_lst[[i]] <- diag(sum(est_names == i))
  }

  for(i in u_names) {
    index_u_name <- which(est_names == i)
    o$value[index_u_name] <- G_lst[[i]] %*% o$value[index_u_name]
  }

  ## Restore the names.
  names(o$value) <- est_names

  ## Combine the matrices for all coefficients, making sure the ordering is correct.
  G <- diag.block.combine(G_lst[u_names])
  ## Delta method to get standard errors.
  o$cov <- G %*% o$cov %*% t(G)
  o$sd <- sqrt(diag(o$cov))

  list(o = o, opt = opt, tmb.output.og = tmb_output_og, esa.partial.derivatives = esa_partial_derivatives)
}

## Builds the fitted object from the prepared fit and the TMB results.
build.fit <- function(spec, tmb){
  out <- outFUN(data.par = spec$data.par,
               data.full = spec$data.full,
               data.trap = spec$data.trap,
               data.traps = spec$data.traps,
               data.mask = spec$data.mask,
               data.dists.thetas = spec$data.dists.thetas,
               detfn = spec$detfn,
               param.og = spec$param.og,
               param.og.4cpp = spec$param.og.4cpp,
               o = tmb$o,
               tmb.output.og = tmb$tmb.output.og,
               opt = tmb$opt,
               name.fixed.par = spec$name.fixed.par,
               name.extend.par = spec$name.extend.par,
               dims = spec$dims,
               DX.full = spec$DX.full,
               DX.mask = spec$DX.mask,
               fix.input = spec$fix.input,
               bucket.info = spec$bucket_info,
               cue.rates = spec$cue.rates,
               mean.cue.rates = spec$mean.cue.rates,
               A = spec$A,
               survey.length = spec$survey.length,
               sound.speed = spec$sound.speed,
               par.extend = spec$par.extend,
               arg.input = spec$arg.input,
               gam.output = spec$gam_output,
               lst.mean.std = spec$lst_mean_std,
               is.scale = spec$is.scale,
               ss.link = spec$ss.link,
               cutoff = spec$cutoff,
               CL = spec$CL,
               esa.partial.derivatives = tmb$esa.partial.derivatives)

  class(out) <- "acre"
  out
}

