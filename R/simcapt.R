#' Simulating SCR data
#'
#' Simulates SCR capture histories and associated additional
#' information in the correct format for use with the function
#' [fit.acre]. If `fit` is provided then no other arguments
#' are required. Otherwise, at least `traps`, `mask`, and
#' `pars` are needed.
#'
#' @param fit an object generated from the model fitting function "fit.acre()" or
#'            the bootstrap process "boot.acre()". If `fit` is provided, then all other
#'             parameters will be taken from the fitted object, so no other parameters need to be provided.
#' @inheritParams read.acre
#' @inheritParams fit.acre
#' @param param A named list. Component names are parameter names, and
#'     each component is the value of the associated parameter. A
#'     value for the parameter `D`, animal density (or call
#'     density, if it an acoustic survey) must always be provided,
#'     along with values for parameters associated with the chosen
#'     detection function and additional information type(s).
#' @param n.sessions An integer for the number of sessions.
#' @param n.rand 
#' @param random.location
#'
#' @export
sim.capt <- function(fit, detfn, param, model = NULL, traps, control.mask = list(), session.cov = NULL, trap.cov = NULL, loc.cov = NULL,
                    dist.cov = NULL, time.loc.cov = NULL, convert.loc2mask = list(), survey.length = NULL, ss.opts = NULL,
                    cue.rates = NULL, n.sessions = NULL, n.rand = 1, random.location = FALSE, sound.speed = 331){
  if(!missing(fit)) {
    ## If 'fit' is provided, get all information from the fitted object.

    detfn <- get.detfn(fit)

    ## Get the data frame that contains the link function for each parameter.
    dat_par <- get.data.param(fit)

    ## Then get the parameter values before back-transforming.
    param <- get.coef(fit)

    ## The "param" argument in the simulation requires back-transformed values.
    param <- param.transform(param, dat_par, back = TRUE)

    n.sessions <- get.dims.tmb(fit)$n.sessions
    par.extend <- get.par.extend(fit)
    traps <- get.trap(fit)
    mask <- get.mask(fit)
    traps <- df.to.list(traps, n.sessions)
    mask <- df.to.list(mask, n.sessions)
    survey.length <- get.survey.length(fit)
    cue.rates <- get.cue.rates(fit)
    sound.speed <- get.sound.speed(fit)
    ss.opts <- get.ss.opts(fit)

  } else {
    stopifnot(all(!missing(detfn), !missing(param), !missing(traps), !is.null(control.mask$buffer)))
    ## Extract n.sessions first, and convert traps into a list if it is not one.
    if(is(traps, 'list')){
      if(!is.null(n.sessions)) warning("'n.sessions' will be ignored as 'traps' is a list.")
      n.sessions <- length(traps)
    } else {
      stopifnot(any(is(traps, 'matrix'), is(traps, 'data.frame')))
      if(is.null(n.sessions)){
        n.sessions <- 1
      }
    }

    traps <- df.to.list(traps, n.sessions)


    ## Create the mask; traps is guaranteed to be a list, so the mask must be a list too.
    control.mask$traps <- traps
    mask <- do.call('create.mask', control.mask)
    ## If a model is assigned, we need to construct "par.extend" for the simulation.
    par.extend <- par.extend.create(loc.cov = loc.cov, mask = mask, dist.cov = dist.cov, time.loc.cov = time.loc.cov,
                                   convert.loc2mask = convert.loc2mask,
                                   session.cov = session.cov, trap.cov = trap.cov)$output
    par.extend$model <- model
  }

  ## write.csv(par.extend$data$mask, 'df_m_sim.csv', row.names = F)
  tem <- sim.data.prepare(detfn, param, par.extend, traps, mask,
                         survey.length, random.location, n.sessions)
  dat_pars <- tem$dat_pars
  dat.density <- tem$dat.density
  dims <- tem$dims
  info.bucket <- tem$info.bucket

  ## For ihd with a categorical covariate, the model does not seem to work well. Here is an example with sim.study('ihd', fit = T):
  ## the true value for D is c(7.3, -0.2, -0.1), with a three-level factor variable.
  ## With the code below, we can confirm that the simulated density on the log scale is really between 7.1 and 7.3, which is correct.
  ## a = attr(mask[[1]], 'area')
  ## log(range(dat.density$D)/a)
  ## However, the estimate of D is just wild; this needs Ben's help.

  if(n.rand == 1){
    tem <- sim.from.param(detfn, dat_pars, dat.density, random.location,
                            dims, info.bucket, ss.opts, cue.rates, sound.speed)

    capt <- tem$o_capt
    o_cue_rates <- tem$o_cue_rates
    animal_locs <- tem$animal_locs

  } else {
    capt <- vector('list', n.rand)
    animal_locs <- vector('list', n.rand)
    o_cue_rates <- vector('list', n.rand)
    pb <- utils::txtProgressBar(1, n.rand, style = 3)
    for(i in 1:n.rand){
      tem <- sim.from.param(detfn, dat_pars, dat.density, random.location,
                                   dims, info.bucket, ss.opts, cue.rates, sound.speed)
      capt[[i]] <- tem$o_capt
      animal_locs[[i]] <- tem$animal_locs
      
      if(!is.null(tem$o_cue_rates)) o_cue_rates[[i]] <- tem$o_cue_rates
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  return(list(capt = capt, sim_cue_rates = o_cue_rates, animal_locs = animal_locs,
              args = list(detfn = detfn, traps = traps, mask = mask, par.extend = par.extend, ss.opts = ss.opts,
                          cue.rates = cue.rates, sound.speed = sound.speed, survey.length = survey.length)))

}

####################################################################################################################################
## Simulate from the assigned parameters.

sim.data.prepare <- function(detfn, param, par.extend, traps, mask, survey.length, random.location, n.sessions){

  if(is.null(survey.length)){
    survey.length <- rep(1, n.sessions)
  } else {
    if(length(survey.length) == 1){
      survey.length <- rep(survey.length, n.sessions)
    } else {
      stopifnot(length(survey.length) == n.sessions)
    }
  }

  stopifnot(detfn %in% c('hn', 'hhn', 'hr', 'th', 'lth', 'ss'))


  ## Fill in basic information.
  n.traps <- sapply(traps, nrow)
  n.masks <- sapply(mask, nrow)
  buffer <- numeric(n.sessions)
  A <- numeric(n.sessions)
  for(s in 1:n.sessions){
    buffer[s] <- attr(mask[[s]], 'buffer')
    A[s] <- attr(mask[[s]], 'area')
  }

  identical_traps <- all(n.traps == n.traps[1])
  identical_masks <- all(n.masks == n.masks[1])


  param.og <- detfn.params(detfn)


  animal.model <- 'mu' %in% names(param)
  is.toa <- 'sigma.toa' %in% names(param)
  is.bearing <- 'kappa' %in% names(param)
  is.dist <- 'alpha' %in% names(param)
  is.ss <- detfn == 'ss'

  param.og <- c(param.og, 'D')
  if(animal.model) param.og <- c(param.og, 'mu')
  if(is.toa) param.og <- c(param.og, 'sigma.toa')
  if(is.bearing) param.og <- c(param.og, 'kappa')
  if(is.dist) param.og <- c(param.og, 'alpha')

  ## Create a session-trap-mask data frame for everything.
  dat_list <- vector('list', n.sessions)
  for(s in 1:n.sessions){

    traps[[s]] <- as.matrix(traps[[s]])
    mask[[s]] <- as.matrix(mask[[s]])

    ## If we do not choose a location randomly inside a mask, then we can
    ## just generate 'dx' and 'theta' here; otherwise, it is not necessary
    ## to generate them here.
    if(!random.location){
      tem_dx <- distances(traps[[s]], mask[[s]])
      tem_theta <- bearings(traps[[s]], mask[[s]])
      dat_list[[s]] <- data.frame(session = rep(s, n.masks[s] * n.traps[s]),
                                 mask = rep(1:n.masks[s], each = n.traps[s]),
                                 trap = rep(1:n.traps[s], n.masks[s]),
                                 dx = as.vector(tem_dx),
                                 theta = as.vector(tem_theta),
                                 mask_size = A[s])
    } else {
      dat_list[[s]] <- data.frame(session = rep(s, n.masks[s] * n.traps[s]),
                                 mask = rep(1:n.masks[s], each = n.traps[s]),
                                 trap = rep(1:n.traps[s], n.masks[s]),
                                 trap_x = rep(traps[[s]][, 1], n.masks[s]),
                                 trap_y = rep(traps[[s]][, 2], n.masks[s]),
                                 x = rep(mask[[s]][,1], each = n.traps[s]),
                                 y = rep(mask[[s]][,2], each = n.traps[s]),
                                 mask_size = A[s])
    }
  }
  dat <- do.call('rbind', dat_list)
  dat <- data.sort(dat, "data.dists.thetas")

  if(is.null(par.extend)){
    name.extend.par <- NULL
    ## If par.extend is not provided, param must be a list of scalars with the "original values" of all parameters.
    for(i in param.og){
      par.value <- param[[i]]
      if(is.null(par.value)) stop(paste0('please provide parameter: ', i))
      dat[[i]] <- par.value
    }
  } else {
    ## If par.extend is provided, param must be a list. If a parameter is extended, then the corresponding element
    ## in "param" should be a vector with all coefficients of its linear predictor; otherwise, the element should be
    ## simply a scalar of the original value.
    if(!all(names(par.extend) %in% c('data', 'model', 'link', 'scale'))) {
      stop("'par.extend' only accepts 'data', 'model', 'link' and 'scale' as input.")
    }

    ## Since at the beginning of the simulation we have no observations at all, we cannot do any scaling.
    ## We recommend that users standardise numeric covariates first.
    if(!is.null(par.extend$scale)){
      warning(paste0("'scale' will be ignored and it will be regarded as FALSE."))
    }

    if(!all(c('data', 'model') %in% names(par.extend))){
      stop("'data' and 'model' must be provided.")
    }

    name.extend.par <- names(par.extend$model)
    if('sigma.b0.ss' %in% name.extend.par) stop('sigma.b0.ss is not supported for parameter extension.')
    if(!all(name.extend.par %in% param.og)){
      stop('one or more parameters assigned in "par.extend" is not valid according to assigned detection function.')
    }


    tem_dat <- dat[,c('session', 'trap', 'mask')]
    input_data <- par.extend$data
    if(!all(names(input_data) %in% c('session', 'trap',
                                     ## 'animal_ID',
                                     'mask'))){
      stop("only 'session', 'trap', or 'mask' level data could be used as input.")
    }
    if(!is.null(input_data$session)){
      stopifnot('session' %in% colnames(input_data$session))
      tem_dat <- merge(tem_dat, input_data$session, by = 'session', all.x = T)
    }
    if(!is.null(input_data$trap)){
      input_data$trap <- extend.dat.check(input_data$trap, 'trap', tem_dat, n.sessions, n.traps, identical_traps)
      tem_dat <- merge(tem_dat, input_data$trap, by = c('session', 'trap'), all.x = T)
    }

    if(!is.null(input_data$mask)){
      input_data$mask <- extend.dat.check(input_data$mask, 'mask', tem_dat, n.sessions, n.masks, identical_masks)
      tem_dat <- merge(tem_dat, input_data$mask, by = c('session', 'mask'), all.x = T)
    }

    ## Include 'x' and 'y' in case they are used.
    data.mask <- as.data.frame(do.call('rbind', mask))
    colnames(data.mask) <- c('x', 'y')
    data.mask$session <- rep(1:n.sessions, n.masks)

    tem_mask_index <- seq(n.masks[1])
    if(n.sessions > 1){
      for(s in 2:n.sessions) tem_mask_index <- c(tem_mask_index, seq(n.masks[s]))
    }
    data.mask$mask <- tem_mask_index

    tem_dat <- merge(tem_dat, data.mask, by = c('session', 'mask'), all.x = T)

    ## Sort it to make sure it has the same order as 'dat'.
    tem_dat <- data.sort(tem_dat, "data.dists.thetas")
    tem_dat$gam.resp <- 1
    for(i in param.og){
      par.value <- param[[i]]
      if(is.null(par.value)) stop(paste0('please provide parameter: ', i))

      if(i %in% name.extend.par){
        foo <- stats::as.formula(paste(c('gam.resp', as.character(par.extend$model[[i]])), collapse = ""))
        tem_model <- mgcv::gam(foo, data = tem_dat, fit = FALSE)
        DX <- tem_model$X
        message(paste0("For extended parameter ", i, ", the corresponding covariates are: ", paste(colnames(DX), collapse = ", "), ".\n",
                   "Please make sure the order of coefficients values in the 'param' for this parameter is correct.\n"))
        if(is.null(par.extend$link[[i]])){
          link <- default.link(i)
        } else {
          link <- par.extend$link[[i]]
        }

        dat[[i]] <- unlink.fun(link, DX %*% par.value)

      } else {
        dat[[i]] <- par.value
      }
    }

  }

  ## D is the animal density, whatever the kind of model.
  dat$D <- dat$D * dat$mask_size
  dims <- list(n.sessions = n.sessions, n.traps = n.traps, n.masks = n.masks, A = A)
  ## Create a data frame with session-mask and 'D' only; we will use this to generate the
  ## simulated numbers of animals.
  dat.density <- dat[, c('session', 'mask', 'x'[random.location], 'y'[random.location], 'D', 'mu'[animal.model])]
  dat.density <- dat.density[!duplicated(dat.density[,c('session', 'mask')]),]
  dat.density <- data.sort(dat.density, 'data.mask')
  dat.density$survey.length <- rep(survey.length, n.masks)
  ## Remove unneeded columns from dat (because they have been transferred to dat.density, or are simply not used).
  dat <- dat[, -which(colnames(dat) %in% c('x', 'y', 'D', 'mu', 'mask_size'))]

  return(list(dat_pars = dat, dat.density = dat.density, dims = dims,
              info.bucket = list(animal.model = animal.model, is.ss = is.ss,
                                is.dist = is.dist, is.bearing = is.bearing,
                                is.toa = is.toa, name.extend.par = name.extend.par)))

}


sim.from.param <- function(detfn, dat.pars, dat.density, random.location, dims, info.bucket,
                          ss.opts, cue.rates, sound.speed, keep.animal.locs=FALSE){
  ## First generate the number of animals for each mask point.
  dat.density$n_animals <- rpois(nrow(dat.density), dat.density$D)
  dat.density <- subset(dat.density, dat.density$n_animals != 0)

  if(nrow(dat.density) == 0){
    ## If there are no simulated animals, return an empty data frame.
    return(data.frame(session = numeric(0), trap = numeric(0), ID = numeric(0)))
  }

  ## Split rows with multiple animals at the same mask point, so we have one row for each individual.
  dat.density <- item.split(dat.density, 'n_animals')

  if(random.location){
      tem_list <- vector('list', dims$n.sessions)
    for(s in 1:dims$n.sessions){
      tem_density <- subset(dat.density, dat.density$session == s)

      ## Calculate the side length of each mask cell. Since the mask points are evenly spread
      ## and each cell is a square, we only need to calculate the average spacing in x or y,
      ## which should be exactly equal to the side length.
      tem <- unique(tem_density$x)
      side_len <- (max(tem) - min(tem))/(length(tem) - 1)

      ## Randomly select a location within the mask cell for each ID.
      tem_density$x <- tem_density$x + runif(nrow(tem_density), -0.5 * side_len, 0.5 * side_len)
      tem_density$y <- tem_density$y + runif(nrow(tem_density), -0.5 * side_len, 0.5 * side_len)

      tem_list[[s]] <- tem_density
    }
    dat.density <- do.call('rbind', tem_list)
  }

  ## Generate random n_calls.
  if(info.bucket$animal.model){
    cue.rates <- NULL
    dat.density$animal_ID <- 1:nrow(dat.density)
    ## Since "n_animals" has been split, each row denotes only one animal, so here we can
    ## simply use mu * survey.length as lambda for the Poisson distribution.
    dat.density$n_calls <- rpois(nrow(dat.density), dat.density$mu * dat.density$survey.length)
  } else {
    ## When cue.rates is not assigned, we regard this as a "cue density model", which means n_calls is equal to n_animals * survey.length.
    if(is.null(cue.rates)){
      ## Since we have split n_animals, each row now represents one animal, so n_animals === 1.
      dat.density$n_calls <- dat.density$survey.length
    } else {
      ## If cue.rates is assigned, we are fitting an "animal density model", so we generate random n_calls based on cue.rates and survey.length.
      cue.rates <- mean(cue.rates)
      dat.density$n_calls <- rpois(nrow(dat.density), cue.rates * dat.density$survey.length)
    }
  }

  ## Split n_calls and assign IDs.
  dat.density <- item.split(dat.density, 'n_calls')
  dat.density$ID <- 1:nrow(dat.density)
  
  if (keep.animal.locs) {
    animal_locs <- dat.density[, c("x", "y")]
  }

  ## Drop unneeded columns.
  dat.density <- dat.density[,-which(colnames(dat.density) %in% c('n_animals', 'n_calls', 'survey.length'))]

  ## Merge dat.density (simulated animals and calls) and dat_pars (values of each parameter) to generate our simulated capture history.
  dat_capt <- merge(dat.density, dat.pars, by = c('session', 'mask'), all.x = T)

  if(random.location){
    ## Re-calculate the distance and theta.
    dat_capt$dx <- sqrt((dat_capt$trap_x - dat_capt$x) ^ 2 +
                            (dat_capt$trap_y - dat_capt$y) ^ 2)

    dat_capt$theta <- bearings.by.vec(dat_capt$trap_x, dat_capt$trap_y,
                                     dat_capt$x, dat_capt$y)
    
    ## Drop unneeded columns.
    dat_capt <- dat_capt[,-which(colnames(dat_capt) %in% c('x', 'y', 'trap_x', 'trap_y'))]
  }

  ## Calculate the probability that each trap detects each call.
  name.param.det <- detfn.params(detfn)
  param.det <- vector('list', length(name.param.det))
  names(param.det) <- name.param.det
  for(i in name.param.det) param.det[[i]] <- dat_capt[[i]]

  if(info.bucket$is.ss){
    ss.link <- ss.opts[['ss.link']]
    if(is.null(ss.link)) ss.link <- 'identity'
    cut_off <- ss.opts[['cutoff']]
    stopifnot(!is.null(cut_off))
    ## The ss model is a little different: if the simulated ss > cutoff, the detection probability is 1,
    ## and 0 otherwise.
    essx <- det.prob(detfn, param.det, dat_capt$dx, ss.link)
    dat_capt[['ss']] <- rnorm(nrow(dat_capt), essx, param.det[["sigma.ss"]])
    dat_capt[['det_prob']] <- ifelse(dat_capt[['ss']] > cut_off, 1, 0)
  } else {
    dat_capt[['det_prob']] <- det.prob(detfn, param.det, dat_capt$dx)
  }

  dat_capt$bincapt <- stats::rbinom(nrow(dat_capt), 1, dat_capt[['det_prob']])
  
  if (keep.animal.locs) {
    animal_locs$bincapt <- 0 
    for (i in 1:nrow(animal_locs)) {
      if (sum(subset(dat_capt, dat_capt$ID == i)$bincapt) > 0) {
        animal_locs$bincapt[i] <- 1
      } 
    }
  }
  
  dat_capt <- subset(dat_capt, dat_capt$bincapt == 1)
  
  

  if(nrow(dat_capt) == 0){
    ## If there are no simulated detections, return an empty data frame.
    return(data.frame(session = numeric(0), trap = numeric(0), ID = numeric(0)))
  }


  output <- data.frame(session = dat_capt$session, ID = dat_capt$ID, occasion = 1, trap = dat_capt$trap)
  if(info.bucket$animal.model) output$animal_ID <- dat_capt$animal_ID
  if(info.bucket$is.ss) output$ss <- dat_capt$ss
  ## Simulate extra info.
  if(info.bucket$is.bearing){
    output[['bearing']] <- 0

    if('kappa' %in% info.bucket$name.extend.par){
      ## If kappa is extended, we need to calculate bearings one by one, since rvm() does not
      ## support vectorised computation with different parameters.
      for(i in 1:nrow(dat_capt)){
        output[i, 'bearing'] <- CircStats::rvm(1, dat_capt[i, 'theta'], dat_capt[i, 'kappa'])
      }
    } else {
      ## If kappa is not extended, we can just take the first value, since they are all the same,
      ## and use vectorised computation.
      k <- dat_capt[['kappa']][1]
      output[['bearing']] <- (dat_capt[['theta']] + CircStats::rvm(nrow(dat_capt), 0, k)) %% (2*pi)
    }
  }

  if(info.bucket$is.dist){
    shape <- dat_capt[['alpha']]
    scale <- dat_capt[['dx']] / shape
    output[['dist']] <- rgamma(nrow(dat_capt), shape = shape, scale = scale)
  }

  if(info.bucket$is.toa){
    output[['toa']] <- rnorm(nrow(dat_capt), dat_capt[['dx']]/sound.speed, dat_capt[['sigma.toa']])
  }

  if(info.bucket$animal.model){
    output <- output[order(output$session, output$animal_ID, output$ID, output$trap), ]
  } else {
    output <- output[order(output$session, output$ID, output$trap), ]
  }

  obj <- list(o_capt = output, o_cue_rates = cue.rates)
  
  if (keep.animal.locs) {
    obj$animal_locs <- animal_locs
  }
  
  return(obj)
}


