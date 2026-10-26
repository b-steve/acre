
get.loc.cov <- function(fit){
  output <- fit$arg_input$loc.cov
  return(output)
}


get.mask <- function(fit){
  output <- fit$args$mask
  return(output)
}

get.mask.from.data <- function(dat){
  output <- dat$mask
  return(output)
}

get.data.trap <- function(fit){
  output <- fit$output.tmb$data.traps
  return(output)
}

get.trap <- function(fit){
  df.traps <- get.data.trap(fit)
  colnames(df.traps)[1:2] <- c('x', 'y')
  n.sessions <- max(df.traps$session)
  output <- vector('list', n.sessions)

  for(s in 1:n.sessions){
    output[[s]] <- as.matrix(subset(df.traps, df.traps$session == s, select = c('x', 'y')))
  }

  return(output)
}

get.trap.from.data <- function(dat){
  output <- dat$traps
  return(output)
}

## Retrieve the capture data (optionally binary) for given animal.id / call.id.
get.capt.by.id <- function(fit, call.id, animal.id=NULL, session=1, return.bincapt=F) {
  if (!(session %in% unique(fit$arg_input$captures$session))) {
    stop("Session ", session, " not found in fit.")
  }

  animal.model <- is.animal.model(fit)
  capt <- fit$arg_input$captures
  capt <- capt[capt$session == session,]

  if (is.null(animal.id) & animal.model) {
    stop("'animal.id' argument must be provided for animal id models")
  }
  if (!animal.model & !is.null(animal.id)) {
    stop("'animal.id' argument should only be provided for animal id models")
  }
  if (!is.null(animal.id) & (length(animal.id) > 1 | !is.numeric(animal.id))) {
    stop("'animal.id' argument should be a numeric vector of length 1")
  }

  ## Make sure call.id is either numeric or exactly "all".
  if (is.null(call.id) ||
      !(is.numeric(call.id) || identical(call.id, "all"))) {
    stop("'call.id' must be a numeric vector or 'all'")
  }
  if (any(duplicated(call.id))) {
    warning("'call.id' argument contains duplicate values, which will be ignored")
    call.id <- unique(call.id)
  }

  if (animal.model) {
    if (!(animal.id %in% capt$animal_ID)) {
      stop(paste("Could not find capture history with animal id:", animal.id, "in session:", session))
    }
    ## Make sure to only grab the appropriate animal's capture history.
    capt <- subset(capt, capt$animal_ID == animal.id)
  }

  if (!all(call.id %in% capt$ID) & !identical(call.id, "all")) {
    missing_calls <- call.id[which(!(call.id %in% capt$ID))]
    missing_calls_msg <- paste(missing_calls, collapse = ", ")
    if (animal.model) {
      stop(paste("Could not find call with 'call.id':", missing_calls_msg,
                 ", for animal with 'animal.id':", animal.id,
                 "in session:", session))
    } else {
      stop(paste("Could not find call with 'call.id':", missing_calls_msg,
                 "in session:", session))
    }
  }

  if (!identical(call.id, "all")) {
    capt <- subset(capt, capt$ID %in% call.id)
  }

  if (return.bincapt) {
    ## Return a data frame with an ID column, and 1 column for each trap.
    n.traps <- nrow(get.trap(fit)[[session]])

    bincapt <- lapply(unique(capt$ID), function(id) {
      ## Grab all rows for this ID.
      sub_capt <- subset(capt, capt$ID == id)
      bincapt  <- numeric(n.traps)
      bincapt[sub_capt$trap] <- 1
      c(id, bincapt)
    })
    bincapt <- do.call(rbind, bincapt)

    ## Convert bincapt from matrix to data frame.
    bincapt <- as.data.frame(bincapt, stringsAsFactors = FALSE)
    colnames(bincapt) <- c("ID", paste0("trap", seq_len(n.traps)))

    return(bincapt)
  } else {
    return(capt)
  }
}

get.par.extend <- function(fit){
  output <- fit$args$par.extend
  return(output)
}

get.par.extend.data <- function(fit){
  output <- get.par.extend(fit)
  output <- output$data
  return(output)
}

get.par.extend.name <- function(fit){
  output <- fit$output.tmb$param.extend
  return(output)
}

get.par.extend.covariate <- function(fit){
  ext_name <- get.par.extend.name(fit)
  ## If no parameter is extended, just return NULL.
  if(is.null(ext_name)) return(ext_name)

  output <- get.par.extend(fit)[['model']]
  for(i in names(output)){
    output[[i]] <- all.vars(output[[i]])
  }
  return(output)
}


get.dims.tmb <- function(fit){
  output <- fit$output.tmb$dims
  return(output)
}

get.cue.rates <- function(fit){
  output <- fit$output.tmb$avg_cue_rates
  return(output)
}

get.survey.length <- function(fit){
  output <- fit$args$survey.length
  return(output)
}

get.infotypes <- function(fit){
  output <- fit$infotypes
  return(output)
}

get.dist.theta <- function(fit){
  output <- fit$output.tmb$data.dists.thetas
  return(output)
}

get.detfn <- function(fit){
  output <- fit$output.tmb$detfn
  return(output)
}


get.data.param <- function(fit){
  output <- fit$output.tmb$param.info.table
  return(output)
}

get.coef <- function(fit){
  output <- fit$output.tmb$coef_link
  return(output)
}

get.param.og <- function(fit){
  output <- fit$output.tmb$param.og
  return(output)
}

get.ss.link <- function(fit){
  output <- fit$output.tmb$ss.link
  return(output)
}


get.ss.opts <- function(fit){
  output <- fit$args$ss.opts
  return(output)
}

get.buffer <- function(fit){
  output <- numeric(fit$n.sessions)
  for(s in 1:fit$n.sessions){
    buffer <- attr(fit$args$mask[[s]], 'buffer')
    if(!is.null(buffer))
      output[s] <- buffer
  }
  return(output)
}

get.sound.speed <- function(fit){
  output <- fit$output.tmb$sound.speed
  return(output)
}


get.gam <- function(fit, param = NULL, which.level = NULL){
  if(!is.null(param)){
    output <- fit$output.tmb$gam_output[[param]]
    if(!is.null(which.level)){
      output <- output[[which.level]]
    }
  } else {
    output <- fit$output.tmb$gam_output
    if(!is.null(which.level)){
      for(i in names(output)){
        output[[i]] <- output[[i]][[which.level]]
      }
    }
  }

  return(output)
}

get.DX.new.gam <- function(mod, newdata){
  var_names <- names(mod$var.summary)
  try({newdata <- newdata[, var_names, drop = FALSE]})
  ## Make sure the new data are in the right data type.
  for(i in var_names){
    if(is.null(newdata[[i]])) stop(paste0('covariate: ', i, " is not provided."))
    if(is(mod$var.summary[[i]], 'factor')) newdata[[i]] <- as.character(newdata[[i]])
    if(is(mod$var.summary[[i]], 'numeric')) newdata[[i]] <- as.numeric(newdata[[i]])
  }

  ## Construct the "old data", making sure it has all levels of all categorical variables.
  olddata <- vector('list', length(var_names))
  names(olddata) <- var_names
  for(i in var_names){
    if(is(mod$var.summary[[i]], 'factor')){
      ## The element var.summary[[i]] is a factor if the explanatory variable is categorical,
      ## and levels() returns a character vector, which is just what we want.
      ## One potential problem: if the gam object was built from data that contains a "factor" class
      ## column in the first place, which has, for example, 5 levels, but in the data this column
      ## only has 4 levels. However, I guess this scenario is extremely unlikely.
      olddata[[i]] <- levels(mod$var.summary[[i]])
    } else if(is(mod$var.summary[[i]], 'numeric')){
      olddata[[i]] <- 1
    } else {
      ## Actually I have never seen any other class; just leave an "in-case" step here.
      stop(paste0('unknown data type in one of gam model objects, the variable name is: ', i,
                  ', please check the relevant data.'))
    }
  }

  ## Adjust the lengths, making sure all elements in this list have the same length.
  len <- max(sapply(olddata, 'length'))

  for(i in var_names) olddata[[i]] <- rep(olddata[[i]], length = len)


  olddata <- as.data.frame(olddata)

  dat <- rbind(olddata, newdata)
  dat$gam.resp <- 1
  output <- mgcv::gam(mod$formula, data = dat, fit = FALSE)$X
  output <- output[(nrow(olddata) + 1):nrow(dat), , drop = FALSE]
  return(output)
}

get.extended.par.value <- function(gam, n.col.full, n.col.mask, par.value.linked, new.covariates,
                                  DX.output = FALSE, matrix.par.value = FALSE){
  if(n.col.full > 1){
    gam.model <- gam$gam_non_mask
    DX_full_new <- get.DX.new.gam(gam.model, new.covariates)
  } else {
    DX_full_new <- matrix(1, ncol = 1, nrow = nrow(new.covariates))
  }

  if(n.col.mask > 0){
    gam.model <- gam$gam_mask
    DX_mask_new <- get.DX.new.gam(gam.model, new.covariates)
    ## Get rid of the first column since the intercept is not here.
    DX_mask_new <- DX_mask_new[, -1, drop = FALSE]
  } else {
    DX_mask_new <- NULL
  }

  DX_new <- cbind(DX_full_new, DX_mask_new)
  if(!matrix.par.value){
    output <- as.vector(DX_new %*% as.matrix(par.value.linked, ncol = 1))
  } else {
    output <- par.value.linked %*% t(DX_new)
  }

  if(DX.output){
    return(list(output = output, DX = DX_new))
  } else {
    return(output)
  }
}

#' Calculate and return coefficients for extended models
#'
#' Returns a list of (nmask) x (ntrap) matrices, for each extended parameter.
#' values correspond to estimated parameter value for either
#' (a) The covariate values used to fit the model or
#' (b) The covariate values provided in newdata
#'
#' Each row corresponds to a trap, and each column a mask point.
#' EXAMPLE: let a model have `sigma` mask-level extended by brand (2
#' traps, 1 "sony" and 1 "panasonic"), and also be mask level
#' extended, say with covariate "noise". Then element `[1, 1]` of the
#' returned `sigma` element, corresponds to the coefficient of `sigma`
#' with "sony" brand, and noise at mask point 1. `[2, 1]` to that of
#' "panasonic" brand, and noise at mask point 1.
#'
#' Note that because the parameter estimates are calculated for the entire mask,
#' over every trap, newdata argument will be recycled (or truncated), to the
#' correct dimensions.
#'
#' @param fit an object generated from the model fitting function
#'    [fit.acre] or the bootstrap process [boot.acre].
#' @param mask a list or a data frame. Generated by [create.mask()].
#' @param traps data frame with numeric columns x and y, corresponding to trap
#'    coordinates.
#' @param newdata data.frame; contains any covariates that will be used
#'    for all extended parameters.
#' @param warn Logical; warn users if they have not passed all covariates in
#'    model specification.
#' @param session an integer value indicating the session of detector array and
#'    individual(s) to be plotted.
#'
#' @return a list of (nmask) x (ntrap) columns, for each extended parameter.
#'          values correspond to estimated parameter value for the given covariate values.
#' @keywords internal
get.par.extend.matrix <- function(fit, mask, traps, session=1, newdata=NULL,
                                  warn = T) {
  ## Make sure the model has an extended parameter.
  par_extended <- !is.null(get.par.extend(fit))
  if (!par_extended) {
    stop("Provided model is not parameter extended")
  }

  ## Make sure the user provides all required covariates in the case of newdata.
  if (!is.null(newdata)) {
    all_covariates <- unlist(get.par.extend.covariate(fit))

    if (!all(all_covariates %in% names(newdata))) {
      if (warn) {
        warning(paste("Missing covariate in 'newdata' argument:",
                   all_covariates[-which(all_covariates %in% names(newdata))]),
                "\nCovariate values from model fitting will be used")
      }
    }
  }

  ext_par_data <- get.par.extend.data(fit)
  n_mask_ext <- nrow(mask)
  n_trap_ext <- nrow(traps)

  ## Select the appropriate session for multi-session covariate data.
  for (level in c("mask", "trap", "session")) {
    if (!is.null(ext_par_data[[level]]) &&
        "session" %in% names(ext_par_data[[level]])) {

      ext_par_data[[level]] <- ext_par_data[[level]][
        ext_par_data[[level]]$session == session,
        ,
        drop = FALSE
      ]
    }
  }

  ext_par_newdata <- NULL

  ## Mask-level covariates.
  if (!is.null(ext_par_data$mask)) {

    if (!is.null(newdata)) {
      for (par in names(newdata)) {
        if (par %in% names(ext_par_data$mask)) {
          ext_par_data$mask[[par]] <- newdata[[par]]
        }
      }
    }

    ext_par_newdata <- ext_par_data$mask[
      rep(seq_len(n_mask_ext), each = n_trap_ext),
      ,
      drop = FALSE
    ]
  }

  ## Trap-level covariates.
  if (!is.null(ext_par_data$trap)) {

    if (!is.null(newdata)) {
      for (par in names(newdata)) {
        if (par %in% names(ext_par_data$trap)) {
          ext_par_data$trap[[par]] <- newdata[[par]]
        }
      }
    }

    trap_ext_data <- ext_par_data$trap[
      rep(seq_len(n_trap_ext),
          length.out = n_mask_ext * n_trap_ext),
      ,
      drop = FALSE
    ]

    if (!is.null(ext_par_newdata)) {
      ext_par_newdata <- cbind(ext_par_newdata, trap_ext_data)
    } else {
      ext_par_newdata <- trap_ext_data
    }
  }

  ## Session-level covariates.
  if (!is.null(ext_par_data$session)) {

    session_covariates <- ext_par_data$session[
      rep(1, n_mask_ext * n_trap_ext),
      ,
      drop = FALSE
    ]

    if (!is.null(ext_par_newdata)) {
      ext_par_newdata <- cbind(ext_par_newdata, session_covariates)
    } else {
      ext_par_newdata <- session_covariates
    }
  }

  ## Calculate the coefficients. Note that this is in long format, so it will have to be
  ## reformatted before returning.
  ext_par_values <- predict(fit, newdata = ext_par_newdata, se.fit = F,
                            confidence = F)

  par_values <- list()
  ext_par_names <- get.par.extend.name(fit)
  ## Reformat the predicted values.
  for (par in ext_par_names) {
    ext_par_values_matrix <- matrix(
      ext_par_values[[par]]$Estimate,
      nrow = nrow(traps),
      ncol = nrow(mask),
      byrow = F
    )
    par_values[[par]] <- ext_par_values_matrix
  }

  return(par_values)
}

## For numeric columns in the data frames in par.extend$data, the mean will be
## used as the default value, and for non-numeric columns, the first element
## will be used as the default value.
## Actually, if the argument fit.acre(par.extend = xxx) sets scale = TRUE,
## then the default values for numeric columns will just be 0; however, if scale = FALSE,
## this function could help.

get.default.covariates <- function(fit){
  dat <- get.par.extend.data(fit)
  dat_name <- names(dat)
  o <- vector('list', length(dat_name))
  names(o) <- dat_name

  for(i in dat_name){
    tem <- dat[[i]]
    tem <- tem[, which(!colnames(tem) %in% c('session','ID','trap','mask')), drop = FALSE]
    tem <- as.data.frame(apply(tem, 2, diy.mean, simplify = FALSE))
    o[[i]] <- tem
  }

  output <- o[[dat_name[1]]]
  if(length(dat_name) > 1){
    for(i in 2:length(dat_name)) output <- cbind(output, o[[dat_name[i]]])
  }
  return(output)
}

get.fixed.par.name <- function(fit){
  output <- fit$output.tmb$param.fix
  if(length(output) == 0) output <- NULL
  return(output)
}

get.sv.for.boot <- function(fit){
  ## Depending on whether any parameter is extended in this model,
  ## we can decide whether to use default covariate values (the mean of each covariate).
  is.extended <- !is.null(get.par.extend.name(fit))
  if(is.extended){
    covaraites_default <- get.default.covariates(fit)
    output <- as.list(coef(fit, 'fitted', new.covariates = covaraites_default))
  } else {
    output <- as.list(coef(fit, 'fitted'))
  }

  df_parm <- get.data.param(fit)
  fixed_par <- names(fit$args$fix)

  ## Avoid extreme start values if the parameter is not fixed.
  for(i in names(output)){
    if(!i %in% fixed_par){
      link <- df_parm[which(df_parm$par == i), 'link']
      if(link == 'logit') output[[i]] <- max(min(0.95, output[[i]]), 0.05)
      if(link == 'log') output[[i]] <- max(0.05, output[[i]])
    }
  }

  return(output)
}

## The simulated detections in a bootstrap, keeping only the auxiliary data used in the fit.
get.capt.for.boot <- function(captures, infotypes, animal.model){
  captures <- captures[, c('session', 'animal_ID'[animal.model], 'ID', 'trap', infotypes), drop = FALSE]
  captures <- data.sort(captures, 'data.full')
  rownames(captures) <- NULL
  captures
}

## Obtain the detections from the acre.data object, with calls (and animals) numbered from 1 within
## each session.
get.capt.for.plot <- function(dat){
  capt <- dat$capt
  animal.model <- "animal_ID" %in% colnames(capt)
  data.capt <- convert.natural.number(capt, animal.model, if(animal.model) "both" else "ID")
  data.capt <- data.sort(data.capt, "data.full")
  return(data.capt)
}

get.boot.res <- function(fit, pars){
  output <- fit$boot$boots
  if(!is.null(pars)){
    name_og <- ori.name(colnames(output))
    output <- output[, which(name_og %in% pars), drop = FALSE]
  }

  return(output)
}


get.boot.res.esa <- function(fit){
  output <- fit$boot$res_esa
  colnames(output) <- gsub("_", "\\.", colnames(output))

  return(output)
}


## Get bias from a single bootstrapped matrix.
get.bias <- function(res, coefs){
  f <- function(x) mean(x, na.rm = T)
  bias <- apply(res, 2, f) - coefs

  return(bias)
}

get.coef.names <- function(pars, pars.extended, data.trap, data.mask) {
  ## If there are no extended parameters, just return the default parameter names.
  if (length(pars.extended) == 0) {
    return(pars)
  }

  data_trap_names <- strsplit(colnames(data.trap), " _ ")
  data_mask_names <- strsplit(colnames(data.mask), " _ ")

  ## Bad code: we should pre-define the length and fill it, rather than grow it with iterations.
  name_output <- c()
  u_pars <- unique(pars)

  for (i in seq_along(u_pars)) {
    if (u_pars[i] %in% pars.extended) {
      index_full <- which(sapply(data_trap_names, function(x) x[1] == u_pars[i]))
      index_mask <- which(sapply(data_mask_names, function(x) x[1] == u_pars[i]))

      for (j in index_full) {
        name_output <- c(name_output, paste0(data_trap_names[[j]], collapse = "."))
      }
      for (j in index_mask) {
        name_output <- c(name_output, paste0(data_mask_names[[j]], collapse = "."))
      }
    }
    else {
      name_output <- c(name_output, u_pars[i])
    }
  }

  return(name_output)
}


###################################################################################################
