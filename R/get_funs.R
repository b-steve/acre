
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

## A simple version of create.capt(); the input is the simulated capture history, if there is any,
## which means nrow(capture) > 0.
get.capt.for.boot <- function(captures, dims, infotypes, animal.model){
  ## Since captures is the simulated data, not the original capture history used in the original model fitting,
  ## dims$n.IDs and dims$animal_ID are not valid here, so do not use them.
  n.sessions <- dims$n.sessions
  n.traps <- dims$n.traps
  ## Deal with irregular IDs and remove unused info types.
  captures <- captures[,c('session', 'animal_ID'[animal.model], 'ID', 'trap', infotypes)]
  ## Convert ID and animal_ID to natural numbers to make it easier to build the head data later.
  captures <- convert.natural.number(captures, animal.model)
  captures$bincapt <- 1

  ## Prepare head data containing all combinations of session-animal_ID[animal.model]-ID-trap.


  tem_data_head <- vector('list', n.sessions)
  for(s in 1:n.sessions){
    capt_session <- subset(captures, captures$session == s)
    if(nrow(capt_session) > 0){
      if(animal.model){
        n.animal <- max(capt_session$animal_ID)
        tem <- vector('list', n.animal)
        for(a in 1:n.animal){
          n.ID <- max(subset(capt_session, capt_session$animal_ID == a)$ID)
          tem[[a]] <- data.frame(session = s, ID = rep(seq(n.ID), each = n.traps[s]),
                                trap = rep(seq(n.traps[s]), n.ID), animal_ID = a)
        }

        tem_data_head[[s]] <- do.call('rbind', tem)

      } else {
        n.ID <- max(capt_session$ID)
        tem_data_head[[s]] <- data.frame(session = s, ID = rep(seq(n.ID), each = n.traps[s]),
                                        trap = rep(seq(n.traps[s]), n.ID))

      }

    } else {
      tem_data_head[[s]] <- data.frame(session = s, ID = NA, trap = 1:n.traps[s])
      if(animal.model) tem_data_head[[s]]$animal_ID <- NA
    }

  }

  data_head <- do.call('rbind', tem_data_head)

  if(animal.model){
    o <- merge(data_head, captures, by = c('session', 'ID', 'trap', 'animal_ID'), all.x = TRUE)
  } else {
    o <- merge(data_head, captures, by = c('session', 'ID', 'trap'), all.x = TRUE)
  }

  o <- data.sort(o, 'data.full')
  for(i in c('bincapt', infotypes)) o[,i] <- ifelse(is.na(o[,i]) & !is.na(o[,'ID']), 0, o[,i])

  if(!animal.model){
    n.IDs <- numeric(n.sessions)
    for(s in 1:n.sessions){
      tem <- subset(captures, captures$session == s)
      if(nrow(tem) > 0) n.IDs[s] <- max(tem$ID)
    }


    ## Prepare the structure of the output.
    if(n.sessions == 1){
      output <- vector('list', length(infotypes) + 1)
      names(output) <- c('bincapt', infotypes)
      for(i in names(output)){
        output[[i]] <- matrix(o[,i], ncol = n.traps, byrow = TRUE)
        rownames(output[[i]]) <- seq(n.IDs)
      }
    } else {
      output <- vector('list', n.sessions)
      for(s in 1:n.sessions){
        output[[s]] <- vector('list', length(infotypes) + 1)
        names(output[[s]]) <- c('bincapt', infotypes)
        tem <- subset(o, o$session == s)
        for(i in names(output[[s]])){
          if(n.IDs[s] > 0){
            output[[s]][[i]] <- matrix(tem[,i], ncol = n.traps[s], byrow = TRUE)
            rownames(output[[s]][[i]]) <- seq(n.IDs[s])
          } else {
            output[[s]][[i]] <- matrix(nrow = 0, ncol = n.traps[s])
          }
        }
      }
    }
  } else {
    output <- o
  }


  return(output)

}

## Obtain the capture history from the acre.data object.
get.capt.for.plot <- function(dat){
  all.types <- c("bearing", "dist", "ss", "toa")

  capt <- dat$capt

  ## According to the current data conversion function, if capt is a data frame, it means
  ## the data are individual-identification-embedded data, i.e. "animal_ID" is included.
  if(is(capt, 'data.frame')){
    ## First convert animal_ID and ID to successive natural numbers.
    data.capt <- convert.natural.number(capt, TRUE, "both")

  } else {
    if("bincapt" %in% names(capt)){
      n.sessions <- 1
      n.traps <- ncol(capt$bincapt)
      n.IDs <- nrow(capt$bincapt)
      extra_info <- all.types[all.types %in% names(capt)]
      is.mrds <- "mrds" %in% names(capt)
    } else {
      n.sessions <- length(capt)
      n.traps <- sapply(capt, function(x) ncol(x$bincapt))
      n.IDs <- sapply(capt, function(x) nrow(x$bincapt))
      extra_info <- all.types[all.types %in% names(capt[[1]])]
      is.mrds <- "mrds" %in% names(capt[[1]])
    }

    tem.data.capt <- vector('list', n.sessions)
    for(i in 1:n.sessions){
      if(n.sessions == 1){
        tem <- capt
      } else {
        tem <- capt[[i]]
      }

      number.row <- n.IDs[i] * n.traps[i]
      tem.df <- data.frame(session = rep(i, number.row), ID = numeric(number.row))
      for(j in c("trap", "bincapt", extra_info)) tem.df[[j]] <- numeric(number.row)
      if(is.mrds) tem.df[,c('mrds_x', 'mrds_y')] <- 0

      if(number.row > 0){
        tem.df$ID <- rep(1:nrow(tem$bincapt), n.traps[i])
        tem.df$trap <- rep(1:n.traps[i], each = n.IDs[i])
        for(k in c('bincapt', extra_info)) tem.df[[k]] <- as.vector(tem[[k]])

        if(is.mrds){
          tem$mrds <- as.data.frame(tem$mrds, stringsAsFactors = FALSE)
          colnames(tem$mrds) <- c('mrds_x', 'mrds_y')
          tem$mrds$ID <- 1:nrow(tem$mrds)
          tem.df <- merge(tem.df, tem$mrds, by = "ID")
        }

      }
      tem.data.capt[[i]] <- tem.df
    }
    data.capt <- do.call("rbind", tem.data.capt)
    data.capt <- convert.natural.number(data.capt, FALSE, "ID")
  }

  data.capt <- data.sort(data.capt, "data.full")
  data.capt <- subset(data.capt, data.capt$bincapt == 1)
  data.capt <- data.capt[, -which(colnames(data.capt) == 'bincapt')]

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

get.coef.names <- function(pars, pars.extended, data.full, data.mask) {
  ## If there are no extended parameters, just return the default parameter names.
  if (length(pars.extended) == 0) {
    return(pars)
  }

  data_full_names <- strsplit(colnames(data.full), " _ ")
  data_mask_names <- strsplit(colnames(data.mask), " _ ")

  ## Bad code: we should pre-define the length and fill it, rather than grow it with iterations.
  name_output <- c()
  u_pars <- unique(pars)

  for (i in seq_along(u_pars)) {
    if (u_pars[i] %in% pars.extended) {
      index_full <- which(sapply(data_full_names, function(x) x[1] == u_pars[i]))
      index_mask <- which(sapply(data_mask_names, function(x) x[1] == u_pars[i]))

      for (j in index_full) {
        name_output <- c(name_output, paste0(data_full_names[[j]], collapse = "."))
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
