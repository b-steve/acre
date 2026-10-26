## Builds the capture data used for fitting: one row per call and detector, with bincapt = 1 for the
## detectors that detected the call, and the auxiliary data at those detectors (zero at the others).
## The detections come from create.capt(), and the numbers of sessions and detectors from the traps.
capture.fun <- function(capt, traps, animal.model){
  all.types <- c("bearing", "dist", "ss", "toa")

  ## If covariate data are included, make sure the data are not all NA.
  ## If we attempt to fit a model where ALL the covariate data are missing
  ## (set to NA), then TMB will provide NaNs for the standard errors.
  for (type in all.types) {
    if (type %in% names(capt)) {
      if (all(is.na(capt[[type]]) | capt[[type]] == 0)) {
        stop(paste0("All the covariate data for \"", type, "\" is NA. Please remove
                  this covariate before converting the data and attempting to fit
                  a model."))
      }
    }
  }

  n.sessions <- length(traps)
  n.traps <- sapply(traps, nrow)
  bucket_info <- all.types[all.types %in% colnames(capt)]
  is.mrds <- animal.model && "mrds_x" %in% colnames(capt)
  ## The animal model keeps only the auxiliary data types that were recorded; the call model has
  ## columns for all of them.
  types <- if(animal.model) bucket_info else all.types
  key.cols <- c("animal_ID"[animal.model], "ID")

  n.IDs <- numeric(n.sessions)
  n.animals <- numeric(n.sessions)
  data.capt <- vector('list', n.sessions)
  for(s in 1:n.sessions){
    det <- capt[capt$session == s, , drop = FALSE]
    if(nrow(det) == 0) next
    ## The calls in this session, in the order of the (sorted) detections.
    call.key <- do.call(paste, c(det[key.cols], sep = "\r"))
    calls <- det[!duplicated(call.key), key.cols, drop = FALSE]
    n.calls <- nrow(calls)
    n.IDs[s] <- n.calls
    if(animal.model) n.animals[s] <- length(unique(calls$animal_ID))
    ## Each call at each detector; the call model numbers its calls 1, 2, ... within the session.
    tem <- data.frame(session = rep(s, n.calls * n.traps[s]))
    if(animal.model){
      tem$animal_ID <- rep(calls$animal_ID, each = n.traps[s])
      tem$ID <- rep(calls$ID, each = n.traps[s])
    } else {
      tem$ID <- rep(seq_len(n.calls), each = n.traps[s])
    }
    tem$trap <- rep(seq_len(n.traps[s]), n.calls)
    tem$bincapt <- 0
    for(j in types) tem[[j]] <- 0
    tem$mrds_x <- 0
    tem$mrds_y <- 0
    row <- (match(call.key, unique(call.key)) - 1) * n.traps[s] + det$trap
    tem$bincapt[row] <- 1
    for(j in bucket_info) tem[[j]][row] <- det[[j]]
    if(is.mrds){
      tem$mrds_x[row] <- det$mrds_x
      tem$mrds_y[row] <- det$mrds_y
    }
    data.capt[[s]] <- tem
  }
  data.capt <- do.call("rbind", data.capt)
  if(is.null(data.capt)){
    data.capt <- data.frame(session = numeric(0), animal_ID = character(0)[animal.model], ID = character(0),
                            trap = numeric(0), bincapt = numeric(0))
  }
  for(j in c(all.types, "mrds_x", "mrds_y")){
    if(!j %in% colnames(data.capt)) data.capt[[j]] <- 0
  }

  if(animal.model){
    data.capt$animal_ID <- as.character(data.capt$animal_ID)
  }
  data.capt$ID <- as.character(data.capt$ID)

  ## The distance model is a gamma distribution, so recorded distances must be positive.
  if("dist" %in% bucket_info){
    detected <- !is.na(data.capt$bincapt) & data.capt$bincapt == 1
    if(any(detected & !is.na(data.capt$dist) & data.capt$dist <= 0)){
      stop("Recorded distances must be positive; use NA for a missing distance.")
    }
  }

  return(list(data.capt = data.capt, dims = list(n.sessions = n.sessions, n.traps = n.traps,
                                                 n.IDs = n.IDs, n.animals = n.animals),
              bucket_info = c(bucket_info, "mrds"[is.mrds])))
}


#################################################################################################

trap.fun <- function(traps, dims){
  stopifnot(any(is(traps, 'list'), is.data.frame(traps), is.matrix(traps)))
  n.sessions <- dims$n.sessions
  n.traps <- dims$n.traps
  ## Convert traps to a list of length n.sessions.
  traps <- df.to.list(traps, n.sessions)
  stopifnot(all(sapply(traps, nrow) == n.traps), all(sapply(traps, ncol) == 2))
  
  
  for(i in 1:n.sessions) {
    traps[[i]] <- as.data.frame(traps[[i]])
    colnames(traps[[i]]) <- c("trap_x", "trap_y")
    traps[[i]]$session <- i
    traps[[i]]$trap <- 1:nrow(traps[[i]])
  }
  return(do.call("rbind", traps))
  
}


###################################################################################

mask.fun <- function(mask, dims, animal.model, data.traps, data.full, bucket.info, local, sound.speed){
  
  stopifnot(any(is(mask, 'list'), is.data.frame(mask), is.matrix(mask)))
  n.sessions <- dims$n.sessions
  n.traps <- dims$n.traps
  n.IDs <- dims$n.IDs
  
  mask <- df.to.list(mask, n.sessions)
  
  n.masks <- sapply(mask, nrow)
  
  A <- numeric(n.sessions)
  buffer <- numeric(n.sessions)
  ## The object dists records the distances between every mask point and every trap;
  ## thetas here records the angle of the line between every mask point and trap.
  ## In the original package it was called "bearings", but there is already
  ## a function with this name, so I changed it.
  dists <- vector(mode = "list", length = n.sessions)
  thetas <- vector("list", length = n.sessions)
  data.dists.thetas <- vector('list', n.sessions)
  data.mask <- vector('list', n.sessions)
  
  for (i in 1:n.sessions){
    traps <- data.traps[data.traps$session == i, c("trap_x", "trap_y")]
    
    if(ncol(mask[[i]]) != 2) {
      stop("each element of 'mask' must be a matrix or data frame with 2 columns")
    }
    
    colnames(mask[[i]]) <- c('x', 'y')
    
    
    data.mask[[i]] <- as.data.frame(mask[[i]], stringsAsFactors = FALSE)
    data.mask[[i]]$session <- i
    data.mask[[i]]$mask <- seq(n.masks[i])
    
    
    mask[[i]] <- as.matrix(mask[[i]])
    dists[[i]] <- distances(as.matrix(traps), mask[[i]])

    thetas[[i]] <- bearings(as.matrix(traps), mask[[i]])
    
    ## A mask point exactly on a detector has a distance of zero, which the distance density
    ## cannot use, and an undefined bearing; as in the location densities, use the smallest
    ## positive distance, and a bearing of zero.
    dx <- as.vector(dists[[i]])
    theta <- as.vector(thetas[[i]])
    theta[dx == 0] <- 0
    dx[dx == 0] <- .Machine$double.eps
    data.dists.thetas[[i]] <- data.frame(session = rep(i, n.masks[i] * n.traps[i]),
                                        mask = rep(1:n.masks[i], each = n.traps[i]),
                                        trap = rep(1:n.traps[i], n.masks[i]),
                                        dx = dx, theta = theta)
    
    
    ## Fill in the area and buffer, if missing.
    if (is.null(attr(mask[[i]], "area"))){
      mask.dists <- distances(mask[[i]], mask[[i]])
      A[i] <- min(mask.dists[mask.dists > 0])^2/10000
    } else {
      A[i] <- attr(mask[[i]], "area")
    }
    
    if (is.null(attr(mask[[i]], "buffer"))){
      buffer[i] <- max(apply(dists[[i]], 2, min))
    } else {
      buffer[i] <- attr(mask[[i]], "buffer")
    }
    attr(mask[[i]], "area") <- A[i]
    attr(mask[[i]], "buffer") <- buffer[i]
  }
  ## The objects data.mask and data.dists.thetas are naturally in the right order here.
  data.dists.thetas <- do.call('rbind', data.dists.thetas)
  data.mask <- do.call('rbind', data.mask)
  
  
  ## Deal with the local argument.
  
  all.which.local <- vector(mode = "list", length = n.sessions)
  all.n.local <- vector(mode = "list", length = n.sessions)
  
  if ("mrds" %in% bucket.info){
    local <- TRUE
    for (i in 1:n.sessions){
      
      if(n.IDs[i] > 0){
        mrds <- data.full[data.full$session == i, c("mrds_x", "mrds_y")]
        mrds <- mrds[!duplicated(mrds),]
        all.which.local[[i]] <- find.nearest.mask(as.matrix(mrds), mask[[i]])
        all.n.local[[i]] <- plyr::laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      }
      
    }
  } else if (local){
    
    for (i in 1:n.sessions){
      if(n.IDs[i] > 0){
        if(animal.model){
          ## For each animal_ID, find the detectors with at least one detection of any call from this animal.

          tem <- data.full[data.full$session == i, c("bincapt", "animal_ID", "trap")]
          tem <- agg.sort(tem, 'bincapt', c('animal_ID', 'trap'), function(x) as.numeric(any(x == 1)))
          bincapt <- tem[, c('animal_ID', 'x')]
          colnames(bincapt) <- c('animal_ID', 'bincapt')
          bincapt <- agg.sort(bincapt, 'bincapt', 'animal_ID', function(x) t(x))
          bincapt <- bincapt[, -1]
        } else {
          bincapt <- data.full[data.full$session == i, c("bincapt","ID")]
          bincapt <- agg.sort(bincapt, 'bincapt', "ID", function(x) t(x))
          bincapt <- bincapt[, -1]
        }
        
        all.which.local[[i]] <- find.local(bincapt, dists[[i]], buffer[i])
        all.n.local[[i]] <- plyr::laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      }
    }
  } else {
    ## If not 'local', these two variables are meaningless, so just assign them something.
    for (i in 1:n.sessions){
      if(n.IDs[i] > 0){
        all.n.local[[i]] <- rep(1, n.IDs[i])
        all.which.local[[i]] <- rep(0, n.IDs[i])
      }
    }
  }

  ## Prepare a data set for ID * mask. It holds the times-of-arrival sums of squares and the local
  ## mask points, so it is only needed with times of arrival or local integration; it has one row
  ## per call and mask point, so it can be large.
  
  data.ID_mask <- vector('list', n.sessions)
  if(local) {
    bucket.info <- c(bucket.info, "local")
    tem.df <- vector('list', n.sessions)
  }
  need.ID_mask <- local || 'toa' %in% bucket.info
  
  ## To avoid confusion, u.id is animal_ID-ID if animal.model, and the "ID" column is
  ## temporarily animal_ID-ID as well, until the later step where we restore it.
  for(i in 1:n.sessions){
    if(n.IDs[i] > 0 && need.ID_mask){
      tem <- subset(data.full, data.full$session == i)

      if(animal.model){
        u.id <- unique(paste(tem[['animal_ID']], tem[['ID']], sep = '---'))
        u.a <- sort(unique(tem$animal_ID))
      } else {
        u.id <- unique(tem[["ID"]])
      }
      
      if(local){
        if(animal.model){
          tem.df[[i]] <- data.frame(session = rep(i, sum(all.n.local[[i]])),
                                   animal_ID = rep(u.a, all.n.local[[i]]),
                                   mask = all.which.local[[i]], local = 1,
                                   stringsAsFactors = F)
        } else {
          tem.df[[i]] <- data.frame(session = rep(i, sum(all.n.local[[i]])),
                                   ID = rep(u.id, all.n.local[[i]]),
                                   mask = all.which.local[[i]], local = 1,
                                   stringsAsFactors = F)
        }
        
        ## Here column 'ID' is "animal_ID---ID" if animal.model; it will be separated later.
        data.ID_mask[[i]] <- data.frame(session = rep(i, n.IDs[i] * n.masks[i]),
                                       ID = rep(u.id, each = n.masks[i]),
                                       mask = rep(1:n.masks[i], n.IDs[i]),
                                       stringsAsFactors = F)
      } else {
        ## Here column 'ID' is "animal_ID---ID" if animal.model; it will be separated later.
        data.ID_mask[[i]] <- data.frame(session = rep(i, n.IDs[i] * n.masks[i]),
                                       ID = rep(u.id, each = n.masks[i]),
                                       mask = rep(1:n.masks[i], n.IDs[i]),
                                       local = 1, stringsAsFactors = F)
      }
      
      ## Assign 'toa_ssq' in this ID_mask data set.
      ## The data.full object must be sorted by 'session', 'ID', 'trap';
      ## it has already been sorted well before this step,
      ## but in case of any disorder, we sort the aggregated result again.
      ## We take care with 'order' because 'aggregate' returns
      ## results based on sorted 'list(ID = xx)' argument.
      if('toa' %in% bucket.info){
        toa <- agg.sort(tem, 'toa', c('animal_ID'[animal.model], 'ID'), function(x) t(x))
        toa.ssq <- make.toa.ssq(toa$x, dists[[i]], sound.speed)
        ## Based on the designed structure of 'data.ID_mask' and the structure
        ## of toa.ssq, it needs to be transposed.
        data.ID_mask[[i]][['toa_ssq']] <- as.vector(t(toa.ssq))
      } else {
        data.ID_mask[[i]][['toa_ssq']] <- 0
      }
    }
  }
  
  data.ID_mask <- do.call('rbind', data.ID_mask)
  
  if(animal.model && need.ID_mask){
    data.ID_mask[['animal_ID']] <- gsub("---(.+)", "", data.ID_mask[['ID']])
    data.ID_mask[['ID']] <- gsub("^(.+)---", "", data.ID_mask[["ID"]])
  }
  
  
  
  if(local){
    tem.df <- do.call('rbind', tem.df)
    if(animal.model){
      data.ID_mask <- merge(data.ID_mask, tem.df, by = c('session', 'animal_ID', 'mask'), all = T)
    } else {
      data.ID_mask <- merge(data.ID_mask, tem.df, by = c('session', 'ID', 'mask'), all = T)
    }
    
    data.ID_mask[['local']] <- ifelse(is.na(data.ID_mask[['local']]), 0, data.ID_mask[['local']])
  }
  
  
  return(list(mask = mask, dists = dists, thetas = thetas, n.masks = n.masks, A = A,
              buffer = buffer, all.which.local = all.which.local, all.n.local= all.n.local,
              data.full = data.full, data.dists.thetas = data.dists.thetas,
              data.mask = data.mask, data.ID_mask = data.ID_mask, bucket_info = bucket.info))
}



###############################################################################

par.extend.fun <- function(par.extend, data.full, data.traps, data.mask, animal.model, dims, fulllist.par, extra.args){
  ## Covariates are at the session, detector or mask level, never at the call level, so the
  ## design columns for the session and detector levels need one row per session and detector.
  data.trap <- data.traps[order(data.traps$session, data.traps$trap), c('session', 'trap')]
  
  n.sessions <- dims$n.sessions
  n.traps <- dims$n.traps
  n.masks <- dims$n.masks
  ## It is possible that the traps/masks remain the same throughout all sessions,
  ## and in this case, users can provide just one data frame of traps/masks.

  identical_traps <- all(n.traps == n.traps[1])
  identical_masks <- all(n.masks == n.masks[1])
  
  ihd_merge <- FALSE
  ## If the user is using the ADMB-version input for the ihd model, there will be an argument named 'ihd.opts'
  ## in extra.args; merge it into 'par.extend'.
  if('ihd.opts' %in% names(extra.args)){
    ## Some basic checks of 'ihd.opts'; this cannot be done in the later check of 'par.extend'
    ## because the user may provide 'par.extend' as well, and inside it there may be
    ## some 'mask'-related extension.
    ihd_data <- extra.args$ihd.opts$covariates
    ihd_scale <- extra.args$ihd.opts$scale
    ihd_model <- extra.args$ihd.opts$model
    stopifnot(!is.null(ihd_data) & !is.null(ihd_model))
    if(is.null(ihd_scale)) ihd_scale <- TRUE
    
    if(is.null(par.extend)){
      par.extend <- list(data = list(mask = ihd_data), model = list(D=ihd_model), scale = ihd_scale)
    } else {
      if(!is.null(par.extend$model$D)){
        stop('"ihd.opts" and Density extension in "par.extend" are both provided, please choose one of them.')
      } else {
        par.extend$model$D <- ihd_model
        par.extend$scale <- ihd_scale
        if(is.null(par.extend$data$mask)){
          par.extend$data$mask <- ihd_data
        } else {
          ihd_data <- covariates.mask.check(ihd_data, n.sessions, n.masks, identical_masks)
          ihd_merge <- TRUE
        }
      }
    }
  }
  
  ## The object fgam is for compatibility with the ADMB model.
  fgam <- NULL
  
  ## Scaling is TRUE by default;
  ## create a list to store the mean and standard deviation of each
  ## component in the formula.
  lst_mean_std <- list()
  if(is.null(par.extend$scale)){
    is.scale <- TRUE
  } else {
    is.scale <- par.extend$scale
  }
  
  namelist.par.extend <- fulllist.par[-which(fulllist.par == 'sigma.b0.ss')]
  
  ## A list recording the number of betas in both data.full and data.mask needed for each parameter.
  param_distribution <- vector('list', length(fulllist.par))
  names(param_distribution) <- fulllist.par
  
  

  if(!is.null(par.extend)){
    ########################################################################################################
    ## Key component begins.
    if(!all(names(par.extend) %in% c('data', 'model', 'link', 'scale'))) {
      stop("'par.extend' only accepts 'data', 'model', 'link' and 'scale' as input.")
    }
    if(!all(c('data', 'model') %in% names(par.extend))){
      stop("'data' and 'model' must be provided.")
    }
    
    name.extend.par <- names(par.extend$model)
    if(!all(name.extend.par %in% namelist.par.extend)){
      stop('one or more parameters indicated in "model" are not supported to be extended currently.')
    }
    
    ## Check the 'data' component first.
    
    ## We could include 'session', 'trap' and 'mask' level covariates,
    ## although some parameters such as 'mu' are in reality only extendable with 'animal_ID'
    ## level covariates; we do not impose any restriction here.
    
    input_data <- par.extend$data
    if(!all(names(input_data) %in% c('session', 'trap', 
                                     ## 'animal_ID',
                                     'mask'))){
      stop("only 'session', 'trap', or 'mask' level data could be used as input.")
    }
    
    data.par.non.mask <- data.frame(session = 1:n.sessions)
    
    ## Check 'session' level data:
    ## a data frame with a 'session' column and session level covariates.
    df.s <- input_data$session
    if(!is.null(df.s)){
      stopifnot(is(df.s, 'data.frame'))
      name.s <- colnames(df.s)
      stopifnot('session' %in% name.s)
      stopifnot(length(name.s) > 1)
      if(any(duplicated(df.s$session))) stop('duplicated "session" input.')
      stopifnot(df.s$session %in% unique(data.full$session))
      var.s <- name.s[-which(name.s == 'session')]
      data.par.non.mask <- merge(data.par.non.mask, df.s, by = 'session', all = TRUE)
    } else {
      var.s <- NULL
    }
    
    ## Check 'trap' level data.
    df.t <- input_data$trap
    if(!is.null(df.t)){
      df.t <- extend.dat.check(df.t, 'trap', data.full, n.sessions, n.traps, identical_traps)
      name.t <- colnames(df.t)

      var.t <- name.t[-which(name.t == 'session'| name.t == 'trap')]
      data.par.non.mask <- merge(data.par.non.mask, df.t, by = 'session', all = TRUE)
    } else {
      var.t <- NULL
    }
    
    
    if(any(c('x', 'y') %in% c(var.s, var.t))){
      stop("'x' and 'y' are reserved column names, please choose another name.")
    }
    
    ## Check 'mask' level data.
    
    ## First check whether any formula in "model" contains "x" or "y".
    is_x <- FALSE
    is_y <- FALSE
    for(i in name.extend.par){
      foo <- stats::as.formula(paste(c('gam.resp', as.character(par.extend$model[[i]])), collapse = ""))
      foo_var <- all.vars(foo[[3]])
      if('x' %in% foo_var) is_x <- TRUE
      if('y' %in% foo_var) is_y <- TRUE
    }
    
    df.m <- input_data$mask
    if(!is.null(df.m)){
      df.m <- extend.dat.check(df.m, 'mask', data.mask, n.sessions, n.masks, identical_masks)
      
      ## When both the old argument 'ihd.opts' and par.extend$data$mask are provided,
      ## we need to merge them after confirming they have no covariates in common.
      if(ihd_merge){
        name1 <- colnames(df.m)
        name1 <- name1[-which(name1 == 'session'| name1 == 'mask')]
        name2 <- colnames(ihd_data)
        name2 <- name2[-which(name2 == 'session'| name2 == 'mask')]
        
        if(any(duplicated(c(name1, name2)))){
          stop('One or more covariates are provided in both "ihd.opts" and "par.extend".')
        }
        
        df.m <- merge(df.m, ihd_data, by = c('session', 'mask'))
      }
      
      name.m <- colnames(df.m)
      var.m <- name.m[-which(name.m == 'session'| name.m == 'mask')]
      data.par.mask <- df.m
      ## If 'x' or 'y' is in a formula but not provided in par.extend$data$mask, then
      ## take values from data.mask.
      if(is_x & !('x' %in% var.m)){
        tem <- data.mask[, c('session', 'mask', 'x')]
        data.par.mask <- merge(data.par.mask, tem, by = c('session', 'mask'))
        var.m <- c(var.m, 'x')
      }
      
      if(is_y & !('y' %in% var.m)){
        tem <- data.mask[, c('session', 'mask', 'y')]
        data.par.mask <- merge(data.par.mask, tem, by = c('session', 'mask'))
        var.m <- c(var.m, 'y')
      }
      
    } else {
      if(is_x | is_y){
        data.par.mask <- data.mask[,c('session', 'mask', 'x'[is_x], 'y'[is_y])]
        var.m <- c('x'[is_x], 'y'[is_y])
      } else {
        var.m <- NULL
        ## The data.par.mask data frame is needed anyway, so create an empty one.
        data.par.mask <- data.frame(session = rep(1:n.sessions, n.masks),
                                   mask = do.call('c', lapply(as.list(n.masks), seq)))
      }
    }
    
    ## Check whether there are duplicated column names across these data sets.
    
    var.ex <- c(var.s, var.t, var.m)
    
    if(any(duplicated(var.ex))) {
      stop('duplicated column name across different level data sets.')
    }
    
    gam_output <- vector('list', length = length(name.extend.par))
    names(gam_output) <- name.extend.par
    ## Check 'model'.
    
    ## Deal with the extended parameters one by one.
    for(i in name.extend.par){
      foo <- stats::as.formula(paste(c('gam.resp', as.character(par.extend$model[[i]])), collapse = ""))
      foo_var <- all.vars(foo[[3]])
      if(length(foo_var) == 0) stop('please input a valid formula except "~ 1"')
      if(!all(foo_var %in% var.ex)) stop(paste0('one or more specifed covariates in "model"',
                                                ' is not provided in "data".'))
      
      ## The parameter sigma_toa is not trap extendable.
      if(i == 'sigma_toa' & any(foo_var %in% var.t)) stop("'sigma_toa' is not extendable in 'trap' level.")
      
      ## In the animal_ID embedded model, the parameter 'mu' is not extendable at the 'trap' level either.
      if(i == 'mu' & any(foo_var %in% var.t)) stop("'mu' is not extendable in 'trap' level.")
      ## In fact there should be many more constraints, e.g. "mu" should not be linked to any 'trap' level covariates.
      ## However, from the code's perspective, there is no need for these constraints, since all parameters
      ## are treated the same as the data structure does not change. However, 'sigma_toa' is different, because 'toa_ssq'
      ## is stored in the data frame "data.ID_matrix", which is not involved in any parameter extension, so
      ## there might be a concern about the workflow later if we enable 'trap' level extension of 'sigma_toa'.
      
      
      tem <- separate.formula(foo, var.m)
      
      ## Here "full" does not mean full information; it does not contain mask information.
      
      foo_full <- tem$foo_full
      foo_mask <- tem$foo_mask
      
      foo_full_var <- all.vars(foo_full[[3]])
      foo_mask_var <- all.vars(foo_mask[[3]])
      
      gam_output[[i]] <- vector('list', 2)
      names(gam_output[[i]]) <- c('gam_non_mask', 'gam_mask')
      
      
      if(length(foo_full_var) > 0){
        ## Avoid any NA.
        tem <- as.data.frame(lapply(data.par.non.mask[, foo_full_var, drop = FALSE], is.na), stringsAsFactors = FALSE)
        index <- apply(tem, 1, any)
        if(sum(index) > 0){
          stop(paste0("NA appears in covariates for parameter ", i, ", a review on the data is recommended."))
        }
        
        
        tem.data <- data.par.non.mask
        tem.data[['gam.resp']] <- 1
        
        
        tem_model <- mgcv::gam(foo_full, data = tem.data, fit = FALSE)
        gam_output[[i]][['gam_non_mask']] <- tem_model
        design.matrix <- tem_model$X
        n_col_full <- ncol(design.matrix)
        design.matrix <- as.data.frame(design.matrix, stringsAsFactors = FALSE)
        colnames(design.matrix) <- paste(i, "_", colnames(design.matrix), sep = " ")
        
        ## Standardise the design matrix if is.scale.
        if(is.scale){
          col_mean <- apply(design.matrix, 2, mean)
          col_sd <- apply(design.matrix, 2, stats::sd)

          lst_mean_std[[i]] <- rbind(col_mean, col_sd)
          nrow_X <- nrow(design.matrix)
          mean_matrix <- matrix(rep(col_mean, nrow_X), nrow = nrow_X, byrow = T)
          ## Set its first column to 0, as we do not scale the intercept column, which should just be 1.
          mean_matrix[,1] <- 0
          sd_matrix <- matrix(rep(col_sd, nrow_X), nrow = nrow_X, byrow = T)
          ## Set its first column to 1, as we do not scale the intercept column, which should just be 1.
          sd_matrix[,1] <- 1
          design.matrix <- (design.matrix - mean_matrix)/sd_matrix
        }
        
        
        design.matrix[['session']] <- tem.data[['session']]
        if(!is.null(tem.data[['trap']])) design.matrix[['trap']] <- tem.data[['trap']]
        if(!is.null(tem.data[['animal_ID']])) design.matrix[['animal_ID']] <- tem.data[['animal_ID']]
        data.par.non.mask <- merge(data.par.non.mask, design.matrix,
                                  by = c('session', 'trap'[!is.null(design.matrix[['trap']])],
                                         'animal_ID'[!is.null(design.matrix[['animal_ID']])]), all = TRUE)
        
      } else {
        ## Whether or not this parameter is extended at other levels, the intercept is always
        ## in this data set. On the other hand, whether or not it is extended at the 'mask' level,
        ## there is no intercept column in data.mask.
        data.par.non.mask[[paste0(i, ' _ (Intercept)')]] <- 1
        n_col_full <- 1
        if(is.scale){
          lst_mean_std[[i]] <- matrix(c(1, 0), nrow = 2, dimnames = list(c('col_mean', 'col_sd'), 
                                                                        paste0(i, ' _ (Intercept)')))
        }
      }
      
      
      if(length(foo_mask_var) > 0){
        tem <- as.data.frame(lapply(data.par.mask[, foo_mask_var, drop = FALSE], is.na), stringsAsFactors = FALSE)
        index <- apply(tem, 1, any)
        if(sum(index) > 0){
          stop(paste0("NA appears in covariates for parameter ", i, ", a review on the data is recommended."))
        }
        
        tem.data <- data.par.mask
        tem.data[['gam.resp']] <- 1
        
        
        ## Delete the first column, which is the intercept.
        tem_model <- mgcv::gam(foo_mask, data = tem.data, fit = FALSE)
        gam_output[[i]][['gam_mask']] <- tem_model
        if(i == 'D') fgam <- tem_model
        design.matrix <- tem_model$X[, -1, drop = FALSE]
        n_col_mask <- ncol(design.matrix)
        design.matrix <- as.data.frame(design.matrix, stringsAsFactors = FALSE)
        colnames(design.matrix) <- paste(i, "_", colnames(design.matrix), sep = " ")
        
        ## Standardise the design matrix if is.scale.
        if(is.scale){
          col_mean <- apply(design.matrix, 2, mean)
          col_sd <- apply(design.matrix, 2, stats::sd)
          
          lst_mean_std[[i]] <- cbind(lst_mean_std[[i]], rbind(col_mean, col_sd))
          
          nrow_X <- nrow(design.matrix)
          mean_matrix <- matrix(rep(col_mean, nrow_X), nrow = nrow_X, byrow = T)
          sd_matrix <- matrix(rep(col_sd, nrow_X), nrow = nrow_X, byrow = T)
          
          ## Here we do not modify the first column because in the mask-related design matrix, the
          ## intercept column has been removed.
          design.matrix <- (design.matrix - mean_matrix)/sd_matrix
        }
        
        design.matrix[['session']] <- tem.data[['session']]
        design.matrix[['mask']] <- tem.data[['mask']]
        data.par.mask <- merge(data.par.mask, design.matrix,
                              by = c('session', 'mask'), all = TRUE)
      } else {
        n_col_mask <- 0
      }
      
      param_distribution[[i]] <- c(n_col_full, n_col_mask)
    }
    
    
    
    
    ## The design matrix generated from gam() has been merged into the data, so the original
    ## columns are redundant.
    var.ex.col.index <- which(colnames(data.par.non.mask) %in% var.ex)
    if(length(var.ex.col.index) > 0){
      data.par.non.mask <- data.par.non.mask[, -var.ex.col.index, drop = FALSE]
    }
    
    var.ex.col.index <- which(colnames(data.par.mask) %in% var.ex)
    if(length(var.ex.col.index) > 0){
      data.par.mask <- data.par.mask[, -var.ex.col.index, drop = FALSE]
    }
    
    
    data.trap <- merge(data.trap, data.par.non.mask, by = c('session', 'trap'[!is.null(data.par.non.mask[['trap']])]),
                       all = TRUE)
    data.mask <- merge(data.mask, data.par.mask, by = c('session', 'mask'))
    
    for(i in fulllist.par){
      if(!(i %in% name.extend.par)){
        data.trap[[paste0(i, ' _ (Intercept)')]] <- 1
        param_distribution[[i]] <- c(1, 0)
      }
    }
    ## Key component ends.
    ########################################################################################################
    
  } else {
    for(i in fulllist.par){
      data.trap[[paste0(i, ' _ (Intercept)')]] <- 1
      param_distribution[[i]] <- c(1, 0)
    }
    name.extend.par <- NULL
    gam_output <- NULL
  }
  
  ## Finally, record the link. Here we only check whether the name is in the full list of parameters,
  ## so technically the user could assign a link function to any parameter here, but this is not recommended.
  df.par <- data.frame(par = character(0), link = character(0), n_col_full = numeric(0),
                      n_col_mask = numeric(0), stringsAsFactors = FALSE)
  
  link <- par.extend$link
  
  if(!is.null(link)){
    if(!all(names(link) %in% fulllist.par)){
      warning("one or more parameters denoted in 'link' are ignored as they are invalid.")
    }
  }
  
  j <- 1
  for(i in fulllist.par){
    df.par[j, 'par'] <- i
    df.par[j, 'n_col_full'] <- param_distribution[[i]][1]
    df.par[j, 'n_col_mask'] <- param_distribution[[i]][2]
    if(i %in% names(link)){
      ## If 'link' is provided, check whether it is a valid input;
      ## currently we only support three link functions.
      stopifnot(link[[j]] %in% c('log', 'logit', 'identity'))
      
      ## Currently we only support the log function as the link function for density
      ## because of the function delta.for.pred() in the file "support_functions.R".
      if(i == 'D') stopifnot(link[[i]] == 'log')
      
      df.par[j, 'link'] <- link[[i]]
    } else {
      ## Set default 'link' if it is not provided.
      df.par[j, 'link'] <- default.link(i)
    }
    j <- j + 1
  }
  
  ## D only accepts log as its link function.
  stopifnot(df.par[which(df.par$par == 'D'), 'link'] == 'log')
  
  data.mask <- data.sort(data.mask, "data.mask")
  
  ## The object fgam is only valid when 'D' is extended at the mask level;
  ## this component is for compatibility with the ADMB model.
  if(is.null(fgam)){
    tem <- data.frame(gam.resp = rep(1, sum(n.masks)))
    fgam <- mgcv::gam(gam.resp~1, data = tem, fit = FALSE)
  }
  
  
  data.trap <- data.trap[order(data.trap$session, data.trap$trap), , drop = FALSE]

  return(list(data.full = data.full, data.trap = data.trap, data.mask = data.mask, data.par = df.par, is.scale = is.scale,
              name.extend.par = name.extend.par, fgam = fgam, gam_output= gam_output,
              lst_mean_std = lst_mean_std))
}

## Construct the Jacobian used to back-transform the scaling and centring done
## on the model covariates.
construct.G.matrix <- function(mean.sd) {
  G <- diag(1 / mean.sd[2,])
  G[1,] <- -1 * mean.sd[1,] / mean.sd[2,]
  ## The intercept needs to be reset manually, as its standard deviation will be
  ## 0, and 1 / standard deviation = infinity.
  G[1,1] <- 1
  
  return(G)
}

## Construct a list of back-transformation matrices that undo any changes
## in scale and shift, where present. Where they are
## not, the identity matrix is returned for the corresponding parameter.
construct.G.matrix.list <- function(par.names, extended.par.names, mean.sd.list, is.scale) {
  u_names <- unique(par.names)
  
  ## Contains the first-derivative matrices for each parameter.
  G_list <- vector('list', length(u_names))
  names(G_list) <- u_names
  
  for(i in u_names) {
    index_u_name <- which(par.names == i)
    np <- length(index_u_name)
    if(i %in% extended.par.names & is.scale) {
      ## Extract the mean_sd matrix. This matrix's first row is the mean, and its second row is the standard deviation.
      ## Column indices correspond to the location of each coefficient within the formula for
      ## that extended parameter.
      mean_sd <- mean.sd.list[[i]]
      
      ## Generate the Jacobian for the given parameter.
      G_list[[i]] <- construct.G.matrix(mean_sd)
    } else {
      ## If the parameter is not extended, there is no need to worry about back-
      ## transforming, as there was never a transformation to begin with.
      G_list[[i]] <- diag(nrow = np)
    }
  }
  
  return(G_list)
}

ss.fun <- function(ss.opts, data.full, data.ID.mask, animal.model, dims, bucket.info, sv, fix){
  cutoff <- ss.opts[["cutoff"]]
  ss.link <- ss.opts[["ss.link"]]

  if("ss" %in% bucket.info){
    if(is.null(ss.opts)) stop("Argument 'ss.opts' is missing.")
    if(is.null(cutoff)) stop("The 'cutoff' component of 'ss.opts' must be specified.")
    ## Warning for unexpected component names.
    if (!all(names(ss.opts) %in% c("cutoff", "ss.link"))){
      warning("Components of 'ss.opts' may only consist of \"cutoff\", \"ss.link\";
            others are being ignored.")
    }
    
    ## Warning for heterogeneity and directional calling models.
    if(any(c("ss.dir", "ss.het") %in% bucket.info)){
      stop("Fitting of models with either heterogeneity in source signal strength
       or directional calling is not yet implemented.")
    }
    
    ## Set the default for ss.link.
    if (is.null(ss.link)){
      ss.link <- "identity"
    } else if (!(ss.link %in% c("identity", "log", "spherical"))){
      stop("Component 'ss.link' in 'ss.opts' must be \"identity\", \"log\", or \"spherical\".")
    }
    

    
    data.ss <- subset(data.full, !is.na(data.full$ID))
    data.no.det.session <- subset(data.full, is.na(data.full$ID))
    
    
    ## Set all capture histories to zero where "ss" < cutoff, and then remove the observations if
    ## there is no capture history for that "ID".
    rem <- data.ss$ss < cutoff
    for(i in c('bincapt', 'ss', 'bearing', 'toa', 'dist')){
      data.ss[[i]] <- ifelse(rem, 0, data.ss[[i]])
    }
    
    keep <- agg.sort(data.ss, 'bincapt', c('session', 'animal_ID'[animal.model], 'ID'), sum)
    
    keep <- subset(keep, keep$x > 0)
    if(nrow(keep) == 0) stop("None of signal received is greater than the cut off threshold.")
    colnames(keep)[which(colnames(keep) == 'x')] <- "keep"
    data.ss <- merge(data.ss, keep, by = c("session", "animal_ID"[animal.model], "ID"), all = TRUE)
    data.ss <- subset(data.ss, !is.na(data.ss$keep))
    
    data.ss <- data.ss[,-which(colnames(data.ss) == "keep")]

    ## Sessions whose detections are all below the cutoff become sessions with no detections,
    ## so they need the same placeholder rows (one per trap, with NA 'ID') as the original ones.
    lost.sessions <- setdiff(unique(data.full$session[!is.na(data.full$ID)]), unique(data.ss$session))
    for(s in lost.sessions){
      tem <- data.full[which(data.full$session == s), ]
      tem <- tem[!duplicated(tem$trap), ]
      for(i in c("animal_ID"[animal.model], "ID", "bincapt", "bearing", "dist", "ss", "toa", "mrds_x", "mrds_y")){
        if(!is.null(tem[[i]])) tem[[i]] <- NA
      }
      data.no.det.session <- rbind(data.no.det.session, tem)
    }

    ## IDs were converted to successive natural numbers before the cutoff was applied; removing
    ## some of them leaves gaps, but TMB uses the ID to index data.full, so convert them again.
    data.ss <- convert.natural.number(data.ss, animal.model, 'both')
    data.ss <- data.sort(data.ss, "data.full")
    data.full <- rbind(data.ss, data.no.det.session)
    data.full <- data.sort(data.full, "data.full")

    ## Refresh "dims"; sessions without any remaining detections should have zero IDs.
    if(animal.model){
      new.n.animals <- agg.sort(data.ss, 'animal_ID', 'session', function(x) length(unique(x)))
      new.n.IDs <- agg.sort(data.ss, c('animal_ID', 'ID'), 'session', function(x) length(unique(x)))
      new.n.animal.call <- agg.sort(data.ss, 'ID', c('session', 'animal_ID'), function(x) length(unique(x)))

      dims$n.animals[] <- 0
      dims$n.animals[new.n.animals$session] <- new.n.animals$x
      dims$n.animal.call <- new.n.animal.call$x
    } else {
      new.n.IDs <- agg.sort(data.ss, 'ID', 'session', function(x) length(unique(x)))
    }
    
    n.removed <- sum(dims$n.IDs) - sum(new.n.IDs$x)
    if(n.removed > 0){
      message(n.removed, " capture history entries have no received signal",
              " strengths above the cutoff and have therefore been removed.\n",
              sep = "")
      ## A logical value recording whether any ID has been removed.
      ID.removed <- TRUE
    } else ID.removed <- FALSE
    
    dims$n.IDs[] <- 0
    dims$n.IDs[new.n.IDs$session] <- new.n.IDs$x

    if(ID.removed && !is.null(data.ID.mask)){
      ## When any ID is deleted due to the ss cutoff, keep the remaining ones (matched on the IDs
      ## before re-conversion), then convert them in the same way as data.full so the two still match.
      if("animal_ID" %in% colnames(data.full)){
        data.ID.mask <- subset(data.ID.mask, paste(data.ID.mask$session, data.ID.mask$animal_ID, data.ID.mask$ID, sep = "-") %in%
                                unique(paste(keep$session, keep$animal_ID, keep$ID, sep = "-")))
      } else {
        data.ID.mask <- subset(data.ID.mask, paste(data.ID.mask$session, data.ID.mask$ID, sep = "-") %in%
                                unique(paste(keep$session, keep$ID, sep = "-")))
      }
      data.ID.mask <- convert.natural.number(data.ID.mask, animal.model, 'both')
      data.ID.mask <- data.sort(data.ID.mask, "data.ID_mask")
    }
    
    ## The modified ss.opts contains much more condensed and regular information.
    ss.opts <- list(ss.link = ss.link,
                   cutoff = cutoff)
    
  } else {
    if (!is.null(ss.opts)){
      warning("Argument 'ss.opts' is being ignored as signal strength is not provided.")
    }
    
    ss.opts <- NULL
    
    if(any(!is.null(sv$b0.ss),
           !is.null(sv$b1.ss),
           !is.null(sv$b2.ss),
           !is.null(sv$sigma.ss),
           !is.null(sv$sigma.b0.ss),
           !is.null(fix$b0.ss),
           !is.null(fix$b1.ss),
           !is.null(fix$b2.ss),
           !is.null(fix$sigma.ss),
           !is.null(fix$sigma.b0.ss))){
      warning(paste0("Not 'ss' model, coresponding parameters specified in",
                     "'sv' or 'fix' will be ignored"))
    }
  }

  return(list(data.full = data.full, dims = dims, data.ID_mask = data.ID.mask,
              ss.opts = ss.opts, bucket_info = bucket.info))
}

############################################################################
CR.SL <- function(cue.rates, survey.length, bucket.info, dims){
  n.sessions <- dims$n.sessions
  
  if (is.null(survey.length)){
    survey.length <- rep(1, n.sessions)
    if (!is.null(cue.rates)){
      warning("The argument `survey.length' is not provided, so it will be assumed to be 1.
              Please ensure `cue.rates' is measured in the same time units.")
    }
  } else {
    if(length(survey.length) != 1){
      if (length(survey.length) != n.sessions){
        stop("The argument `survey.length' must have a value for each session.")
      }
    } else {
      survey.length <- rep(survey.length, n.sessions)
    }

  }
  
  if (!is.null(cue.rates)){
    bucket.info <- c(bucket.info, "freqs")
    mean.cue.rates <- mean(cue.rates)
  } else {
    mean.cue.rates <- 1
  }
  
  return(list(bucket_info = bucket.info, mean.cue.rates = mean.cue.rates, survey.length = survey.length))
  
}


############################################################################

param.detfn.fun <- function(animal.model, sv, fix, bounds, name.extend.par, detfn, data.full, data.trap, data.mask,
                           data.par, ss.opts, bucket.info, fulllist.par, A, buffer,
                           survey.length, dims, sv.link){
  ## Check bounds, sv, fix.
  if(!is.null(bounds)){
    if(!is(bounds, 'list')) stop ('"bounds" must be a list.')
    if(!all(names(bounds) %in% fulllist.par)) stop('one or more element in "bounds" is not valid.')
    if(!all(sapply(bounds, function(x) nrow(as.matrix(x))) == 2)){
      stop('each element of "bounds" should be a vector with two values, or a matrix with two rows.')
    }
  } else {
    bounds <- vector('list', 0)
  }
  
  if(!is.null(sv)){
    if(!is(sv, 'list')) stop('"sv" must be a list.')
    if(!all(names(sv) %in% fulllist.par)) stop('one or more element in "sv" is not valid.')
  } else {
    sv <- vector('list', 0)
  }
  
  if(!is.null(fix)){
    if(!is(fix, 'list')) stop('"fix" must be a list.')
    if(!all(names(fix) %in% fulllist.par)) stop('one or more element in "fix" is not valid.')
  } else {
    fix <- vector('list', 0)
  }
  
  
  
  ## The detection function is the basis for choosing parameters, so settle this first.
  
  if(is.null(detfn)){
    ## The detfn for 'ss' is a little different from the original, and may need
    ## to be restored to the original at the very end of the 'fit' function.
    if("ss" %in% bucket.info){
      detfn <- 'ss'
    } else {
      detfn <- "hn"
    }
  } else {
    if("ss" %in% bucket.info){
      if(detfn != 'ss'){
        warning('signal strenth is provided, therefore the specified "detfn" is ignored')
      }
      detfn <- 'ss'

      
    } else {
      if(detfn == 'ss'){
        warning('signal strenth is not provided, therefore the specified "detfn" is replaced to "hn" (half normal),')
        detfn <- 'hn'
      } else {
        if(!detfn %in% c('hn', 'hhn', 'hr', 'th', 'lth')){
          warning('invalid "detfn" is specified, the default will be used instead')
          detfn <- 'hn'
        }
      }
    }
  }
  
  ## Then decide the original names of the parameters that will be input into TMB.
  param.og <- detfn.params(detfn)

  
  if('bearing' %in% bucket.info){
    param.og <- c(param.og, 'kappa')
  }
  
  if('dist' %in% bucket.info){
    param.og <- c(param.og, 'alpha')
  }
  
  if('toa' %in% bucket.info){
    param.og <- c(param.og, 'sigma.toa')
  }
  
  param.og <- c(param.og, "D")
  
  if(animal.model) param.og <- c(param.og, "mu")
  
  ## Re-order param.og to make sure its order is consistent with "fulllist.par.generator()".
  param.og <- fulllist.par.generator()[fulllist.par.generator() %in% param.og]
  
  ## Generate the DATA_MATRIX for every parameter,
  ## and in the loop, also do some basic conflict checking.
  data.full <- data.sort(data.full, "data.full")
  
  
  design.matrices.full <- vector('list', length(fulllist.par))
  names(design.matrices.full) <- fulllist.par
  
  tem <- strsplit(colnames(data.trap), " _ ")
  tem <- lapply(tem, function(x) x[[1]])
  clean.names.full <- do.call('c', tem)
  
  design.matrices.mask <- vector('list', length(fulllist.par))
  names(design.matrices.mask) <- fulllist.par
  
  tem <- strsplit(colnames(data.mask), " _ ")
  tem <- lapply(tem, function(x) x[[1]])
  clean.names.mask <- do.call('c', tem)
  
  for(i in fulllist.par){
    index_par <- which(data.par$par == i)
    
    if(data.par[index_par, "n_col_mask"] > 0){
      design.matrices.mask[[i]] <- as.matrix(data.mask[, which(clean.names.mask == i), drop = FALSE])
    } else {
      design.matrices.mask[[i]] <- matrix(0, ncol = 2, nrow = 2)
    }
    
    ## Because the intercept is in "data.trap", this design matrix always has something in it.
    design.matrices.full[[i]] <- as.matrix(data.trap[, which(clean.names.full == i), drop = FALSE])
    
    if(i %in% param.og){
      ## Deal with "fix" first.
      if(i %in% name.extend.par){
        ## A valid extended parameter should not be fixed.
        if(!is.null(fix[[i]])){
          warning(paste0("parameter ", i, " is going to be extended, therefore, it cannot be fixed."))
          fix[[i]] <- NULL
        }
      } else {
        ## If this valid parameter is not extended, check whether its fixed value is in range.
        if(!is.null(fix[[i]])){
          tem <- param.range.validate(i, fix[[i]])
          if(!tem) stop(paste0("fixed value for parameter ", i, " is out of valid range."))
        }
      }
      ## Then deal with "sv" and "bounds".
      ## If a parameter is extended, the start value will be its intercept beta's start value,
      ## and the bounds of an extended parameter will be assigned to its intercept as well.
      if(!is.null(bounds[[i]])){
        ## The lower and upper bounds should both be in the valid range of this parameter.
        tem1 <- param.range.validate(i, ifelse(bounds[[i]][1] == 0, 1e-15, bounds[[i]][1]))
        tem2 <- param.range.validate(i, ifelse(bounds[[i]][2] == 0, 1e-15, bounds[[i]][2]))
        if(!all(tem1, tem2)) stop(paste0("bounds for parameter ", i, " is out of valid range."))
        if(bounds[[i]][2] <= bounds[[i]][1]) stop(paste0("invalid range for parameter ", i))
        
        if(!is.null(sv[[i]])){
          tem <- param.range.validate(i, sv[[i]])
          if(!tem) stop(paste0("start value for parameter ", i, " is out of valid range."))
          
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            stop(paste0("the start value of ", i, " is out of the specified corresponding bounds."))
          }
        } else {
          ## If sv is not assigned, assign a default value.
          sv[[i]] <- default.sv(i, info = list(data.full = data.full, data.mask = data.mask,
                                              buffer = buffer, A = A, ss.opts = ss.opts,
                                              detfn = detfn, dims = dims, animal.model = animal.model,
                                              survey.length = survey.length, param.og = param.og,
                                              sv = sv))
          
          ## Since bounds are provided, we need the default start value to lie within them;
          ## if the default start value is outside the specified bounds, re-assign the start value.
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            sv[[i]] <- 0.5 * (bounds[[i]][1] + bounds[[i]][2])
          }
        }
        
        ## Check whether the fixed value is inside the assigned bounds.
        ## Since we have already checked the validity of the fixed value, there is no need to check it here as for sv;
        ## if bound is NULL, there is no need to check the fixed value against the default bounds.
        if(!is.null(fix[[i]])){
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            stop(paste0("the start value of ", i, " is out of the specified corresponding bounds."))
          }
        }
        
      } else {
        bounds[[i]] <- default.bounds(i)
        
        if(!is.null(sv[[i]])){
          tem <- param.range.validate(i, sv[[i]])
          if(!tem) stop(paste0("start value for parameter ", i, " is out of valid range."))
          ## I'm not sure what we should do if the assigned start value is outside the default bounds;
          ## discuss with Ben later.
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            warning(paste0("specified start value of '", i, "' is extreme, a carefull checking ",
                           "is recommended."))
            bounds[[i]][1] <- -2 * abs(sv[[i]])
            bounds[[i]][2] <- 2 * abs(sv[[i]])
          }
          
        } else {
          ## If sv is not assigned, assign a default value.
          
          sv[[i]] <- default.sv(i, info = list(data.full = data.full, data.mask = data.mask,
                                              buffer = buffer, A = A, ss.opts = ss.opts,
                                              detfn = detfn, dims = dims, animal.model = animal.model,
                                              survey.length = survey.length, param.og = param.og,
                                              sv = sv))
          
        }
      }
    } else {
      ## A start value is needed for every parameter, regardless of whether it is used.
      if(!is.null(sv[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'sv' will be ignored because ",
                       "the detect function is incompatible"))
      }
      sv[[i]] <- 0.5
      
      ## If this parameter is not used in the model but is in "fix", remove it from "fix"
      ## by assigning NULL to it.
      if(!is.null(fix[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'fix' will be ignored because ",
                       "the detect function is incompatible"))
        fix[[i]] <- NULL
      }
      
      if(!is.null(bounds[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'bounds' will be ignored because ",
                       "the detect function is incompatible"))
      }
      
      bounds[[i]] <- c(0.3, 0.7)
      
      ## If this unused parameter is specified to be extended, everything relating to the extension will be ignored.
      if(i %in% name.extend.par){
        warning(paste0("parameter ", i, " specified to be extended, however, this parameter is not used ",
                       "in the model, everything about extension will be ignored"))
        ## Correct all of the relevant information.
        design.matrices.full[[i]] <- design.matrices.full[[i]][, 1, drop = FALSE]
        design.matrices.mask[[i]] <- matrix(0, ncol = 2, nrow = 2)
        data.par[index_par, 'n_col_mask'] <- 0
        data.par[index_par, 'n_col_full'] <- 1
        name.extend.par <- name.extend.par[-which(name.extend.par == i)]
      }
      
      
    }
  }
  
  ## Construct the parameter input, which must be transformed;
  ## the elements do not need to be named, especially the betas for extended parameters.
  
  ## The sv.input list contains vectors as elements;
  ## in each vector, the first 'n_col_full' numbers are for the relevant columns in
  ## data.full, and the rest for the relevant columns in data.mask.
  ## The fix.input list contains scalars as elements, because an extended parameter
  ## cannot be fixed.
  sv.input <- vector('list', length(fulllist.par))
  names(sv.input) <- fulllist.par
  fix.input <- vector('list', length(fulllist.par))
  names(fix.input) <- fulllist.par
  ## The bounds.input list contains matrices as elements.
  bounds.input <- vector('list', length(fulllist.par))
  names(bounds.input) <- fulllist.par
  
  
  for(i in fulllist.par){
    index_par <- which(data.par$par == i)
    length_param <- data.par[index_par, 'n_col_full'] + data.par[index_par, 'n_col_mask']
    link <- data.par[index_par, 'link']
    
    
    ## There is a small problem that needs attention: if the user specifies a parameter to be fixed at zero
    ## and its link function is 'log', then the real input of the fixed value in TMB will be log(1e-20),
    ## which is not exactly zero. But this problem can be solved by rounding the output to 7 digits,
    ## so in the output this parameter will be zero again.
    if(!is.null(fix[[i]])) {
      fix.input[[i]] <- link.fun(link = link, value = fix[[i]])
      ## If a parameter is fixed, we set its start value to its fixed value;
      ## since only non-extended parameters can be fixed, we don't need to
      ## worry about its length, as it is always a scalar.
      sv.input[[i]] <- fix.input[[i]]
    } else {
      sv.input[[i]] <- numeric(length_param)
      ## Since numeric(n) defaults to 0, we only need to assign the
      ## first element, whether or not it is extended.
      sv.input[[i]][1] <- link.fun(link = link, value = sv[[i]])
    }
    
    ####################################################################################
    ## This is mostly for development purposes; its use is not recommended.
    if(!is.null(sv.link[[i]])) sv.input[[i]] <- sv.link[[i]]
    ####################################################################################
    
    ## Set the default bounds of the extended betas as the starting point of our bounds matrix,
    ## and replace the first one (or the only one) with the transformed assigned bounds.
    if(link == 'identity'){
      bounds.input[[i]] <- matrix(rep(c(-1e8, 1e8), length_param),
                                 nrow = 2)
    } else {
      bounds.input[[i]] <- matrix(rep(c(-log(1e20), log(1e20)), length_param),
                                 nrow = 2)
    }
    
    bounds.input[[i]][1, 1] <- link.fun(link = link, value = bounds[[i]][1])
    bounds.input[[i]][2, 1] <- link.fun(link = link, value = bounds[[i]][2])
  }

  
  

  
  
  return(list(design.matrices.full = design.matrices.full,
              design.matrices.mask = design.matrices.mask,
              name.extend.par = name.extend.par,
              data.par = data.par,
              sv.input = sv.input,
              fix.input = fix.input,
              bounds.input = bounds.input,
              detfn = detfn,
              param.og = param.og))
  
}



gr.free.o.restore <- function(fn, opt, H, parameters, param.og.4cpp, n.sessions){
  parameters <- parameters[param.og.4cpp]
  name_fitted_par <- unique(names(opt$par))
  for(i in name_fitted_par) parameters[[i]] <- val(opt$par[which(names(opt$par) == i)])
  o_val <- list.2vector.4value(parameters)
  esa <- numeric(n.sessions)
  names(esa) <- rep('esa', n.sessions)
  o_val <- c(o_val, esa)
  
  len <- length(o_val)
  o_val_names <- names(o_val)
  i_fitted_par <- which(o_val_names %in% name_fitted_par)
  
  
  sd_fixed <- sqrt(diag(H))
  o_sd <- numeric(len)
  o_sd[i_fitted_par] <- sd_fixed
  
  o_cov <- matrix(0, nrow = len, ncol = len)
  o_cov[i_fitted_par, i_fitted_par] <- H
  
  o <- vector('list', 5)
  names(o) <- c('value', 'sd', 'cov', 'cov.fixed', 'gradient.fixed')
  o$value <- o_val
  o$sd <- o_sd
  o$cov <- o_cov
  o$cov.fixed <- H
  o$gradient.fixed <- matrix(numDeriv::grad(fn, opt$par), nrow = 1)
  return(o)
}

outFUN <- function(data.par, data.full, data.trap, data.traps, data.mask, data.dists.thetas, detfn, param.og, param.og.4cpp, o, tmb.output.og, opt,
                  name.fixed.par, name.extend.par, dims, DX.full, DX.mask, fix.input, bucket.info, cue.rates, mean.cue.rates, A,
                  survey.length, sound.speed, par.extend, arg.input, fgam, gam.output, lst.mean.std, is.scale, ss.link, cutoff,
                  CL, esa.partial.derivatives){
  ###################################################################################################################
  ## Sort out the output for the function.
  out <- vector('list', 36)
  names(out) <- c("fn", "coefficients", "coeflist", "se", "loglik", "opt_output", "maxgrad", "cor", "vcov", "npar", "npar_re",
                 "npar_sdrpt", "npar_rep", "npar_total", "hes", "eratio", "D.mask", "mm.ihd", "args", "n.sessions",
                 "fit.types", "infotypes", "detpars", "suppars", "D.betapars", "phases", "par.links", "par.unlinks",
                 "fit.ihd", "re.detfn", "fit.freqs", "first.calls", "model.formula", "fgam", "all.covariates",
                 "output.tmb")
  ## Create the output for the TMB model.
  out[['output.tmb']] <- vector('list', 22)
  names(out[['output.tmb']]) <- c('coef_link', 'se_link', 'DX_standardized', 'detfn', 'param.og', 'param.extend', 'param.fix',
                                 'param.info.table', 'data.traps', 'data.full', 'data.mask', 'data.dists.thetas', 'dims',
                                 'avg_cue_rates', 'sound.speed', 'area_unit', 'survey.length', 'ss.link', 'cutoff', 'gam_output',
                                 'tmb_output_og', 'mean_std_for_standardizing')
  
  ## Give each parameter an index, to make it easier to find in "data.par".
  par.id <- 1:nrow(data.par)
  names(par.id) <- data.par$par
  
  param.with.esa <- c(param.og, 'esa')
  param.og.4cpp <- c(param.og.4cpp, 'esa')
  o_value <- vector('list', length(param.with.esa))
  names(o_value) <- param.with.esa
  o_sd <- vector('list', length(param.with.esa))
  names(o_sd) <- param.with.esa
  for(i in 1:length(param.with.esa)){
    o_value[[param.with.esa[i]]] <- o$value[which(names(o$value) == param.og.4cpp[i])]
    o_sd[[param.with.esa[i]]] <- o$sd[which(names(o$value) == param.og.4cpp[i])]
  }
  
  
  ############################################################################################################
  ## The 1st component: "fn".
  out$fn <- 'secr_tmb'
  
  
  #############################################################################################################
  ## The 2nd component: "coefficients"; "esa.[session]" is not included yet.
  ## This can be generated from "name.fited.pars", "name.fixed.par" and "name.extend.par".
  
  ## Create a vector of parameter names used in various output components,
  ## to match the output of the current ADMB model.
  
  ## Here "name.fited.pars" is used instead of "param.with.esa" to make sure the order
  ## of parameters is the same as in the TMB model, and "esa" is not in "param.og";
  ## if nothing goes wrong, "name.fited.pars" should be exactly the same as "param.with.esa".
  name.fited.pars <- unique(names(o$value))
  name.fited.pars <- gsub("_", "\\.", name.fited.pars)
  name_output <- vector('list', length(name.fited.pars))
  names(name_output) <- name.fited.pars
  
  ## Names for the TMB output; "esa" is kept entirely separate from the other parameters.
  name_output_tmb <- name_output
  name_output_tmb$esa <- NULL
  
  name_split_full <- strsplit(colnames(data.trap), " _ ")
  name_split_mask <- strsplit(colnames(data.mask), " _ ")
  
  ## The object out_coef is temporary, and will be assigned to out$coefficients later.
  out_coef <- vector('list', length(name.fited.pars))
  names(out_coef) <- name.fited.pars
  ## Take 'esa' out of the TMB output.
  coef_link_tmb <- out_coef
  coef_link_tmb$esa <- NULL
  
  for(i in name.fited.pars){
    ## The 'esa' is not in the original parameter list, so deal with it separately.
    if(i != 'esa'){
      ## First, sort out the names for the output coefficients.
      index_full <- which(sapply(name_split_full, function(x) x[1] == i))
      index_mask <- which(sapply(name_split_mask, function(x) x[1] == i))
      name_output[[i]] <- c(colnames(data.trap)[index_full], colnames(data.mask)[index_mask])
      name_output_tmb[[i]] <- name_output[[i]]
      name_output[[i]] <- gsub(" _ ", ".", name_output[[i]])
      ## If i is "D", then keep the .(Intercept) whatever the model.
      if(i != "D" & !i %in% name.extend.par){
        name_output[[i]] <- gsub("\\.\\(Intercept\\)", "", name_output[[i]])
      }
      
      if(!i %in% name.fixed.par){
        out_coef[[i]] <- o_value[[i]]
        
        if(is.null(out_coef[[i]])){
          out_coef[[i]] <- fix.input[[i]]
        }
        
        names(out_coef[[i]]) <- paste(name_output[[i]], 'link', sep = "_")
        coef_link_tmb[[i]] <- out_coef[[i]]
        names(coef_link_tmb[[i]]) <- name_output_tmb[[i]]
        
        if(i == "D"){
          tem <- out_coef[[i]]
        } else {
          if(!i %in% name.extend.par){
            tem <- unlink.fun(data.par[par.id[i], 'link'], out_coef[[i]])
          } else {
            tem <- out_coef[[i]]
          }
        }
        
        names(tem) <- name_output[[i]]
        out_coef[[i]] <- c(out_coef[[i]], tem)
        
        if(i == 'D' & !i %in% name.extend.par){
          tem <- unlink.fun(data.par[par.id[i], 'link'], tem)
          names(tem) <- 'D'
          out_coef[[i]] <- c(out_coef[[i]], tem)
        }
      } else {
        coef_link_tmb[[i]] <- fix.input[[i]]
        names(coef_link_tmb[[i]]) <- name_output_tmb[[i]]
      }
    } else if (i == 'esa') {
      name_output[[i]] <- paste('esa', 1:dims$n.sessions, sep = ".")
      out_coef[[i]] <- o_value[[i]]
      names(out_coef[[i]]) <- name_output[[i]]
    }
  }
  
  out$coefficients <- do.call('c', out_coef)
  names(out$coefficients) <- do.call('c', lapply(out_coef, function(x) names(x)))
  out[['output.tmb']][['coef_link']] <- coef_link_tmb
  #######################################################################################################
  ## The 3rd component: "coeflist".
  
  out$coeflist <- vector('list', 0)
  
  D_pars <- name_output$D
  det_pars <- name.fited.pars[which(!name.fited.pars %in% c('D', 'kappa', 'alpha', 'sigma.toa', 'esa'))]
  sup_pars <- name.fited.pars[which(name.fited.pars %in% c('kappa', 'alpha', 'sigma.toa'))]
  
  tem <- val(o_value[['D']])
  for(i in 1:length(name_output$D)){
    tem_name <- paste0('D_betapars_link[', i, ']')
    out$coeflist[[tem_name]] <- tem[i]
  }
  
  for(i in 1:length(det_pars)){
    if(!det_pars[i] %in% name.fixed.par){
      tem <- val(o_value[[det_pars[i]]])
      tem_name <- paste0('detpars_link[', i, ']')
      out$coeflist[[tem_name]] <- tem
    }
  }
  
  if(length(sup_pars) > 0){
    for(i in 1:length(sup_pars)){
      if(!sup_pars[i] %in% name.fixed.par){
        tem <- val(o_value[[sup_pars[i]]])
        tem_name <- paste0('suppars_link[', i, ']')
        out$coeflist[[tem_name]] <- tem
      }
    }
  }
  
  
  for(i in 1:dims$n.sessions){
    tem <- val(o_value[['esa']][i])
    tem_name <- paste0('esa[', i, ']')
    out$coeflist[[tem_name]] <- tem
  }
  
  ## Deal with "D" and "D.mask[s]".
  if(!'D' %in% name.extend.par){
    out$coeflist[['D']] <- unlink.fun(data.par[par.id['D'], 'link'], val(o_value[['D']]))
    for(i in 1:dims$n.sessions){
      tem_name <- paste0('D_mask[', i, ']')
      out$coeflist[[tem_name]] <- rep(out$coeflist[['D']], dims$n.masks[i])
    }
  } else {
    tem_n_col <- data.par[par.id['D'], c('n_col_full', 'n_col_mask')]
    D_est_tmb_og <- tmb.output.og$est[which(names(tmb.output.og$est) == 'D')]
    tem_full <- DX.full[['D']] %*% val(D_est_tmb_og)[1:tem_n_col[1, 1]]
    if(tem_n_col[1, 2] > 0){
      tem_mask <- DX.mask[['D']] %*% val(D_est_tmb_og)[(tem_n_col[1, 1] + 1):sum(tem_n_col)]
    } else {
      tem_mask <- rep(0, sum(dims$n.masks))
    }
    
    for(i in 1:dims$n.sessions){
      tem_name <- paste0('D_mask[', i, ']')
      ## D is only extendable at the session level, not the trap level, so we can take the first row of each session
      ## in 'data.trap', which is tem_full below.
      out$coeflist[[tem_name]] <- tem_full[which(data.trap$session == i)[1]] +
        tem_mask[which(data.mask$session == i)]
      out$coeflist[[tem_name]] <- unlink.fun(data.par[par.id['D'], 'link'],
                                            out$coeflist[[tem_name]])
    }
  }
  
  
  
  ####################################################################################################
  ## The 4th component: "se".
  
  out_se <- vector('list', length(name.fited.pars))
  names(out_se) <- name.fited.pars
  
  ## The TMB output does not include 'esa'.
  se_tmb <- out_se
  se_tmb$esa <- NULL
  
  for(i in name.fited.pars){
    if(i != 'esa'){
      if(!i %in% name.fixed.par){
        ## Then input the estimated values.
        out_se[[i]] <- o_sd[[i]]
        names(out_se[[i]]) <- paste(name_output[[i]], 'link', sep = "_")
        
        se_tmb[[i]] <- o_sd[[i]]
        names(se_tmb[[i]]) <- name_output_tmb[[i]]
        if(i == "D"){
          ## For "D", treat it as extended whether or not it actually is.
          tem <- out_se[[i]]
        } else {
          if(!i %in% name.extend.par){
            tem <- delta.fun(data.par[par.id[i], 'link'], out_se[[i]], o_value[[i]])
          } else {
            tem <- out_se[[i]]
          }
        }
        names(tem) <- name_output[[i]]
        out_se[[i]] <- c(out_se[[i]], tem)
        
        
        ## For homogeneous "D" models, there is an element in the vector named "D".
        if(i == 'D' & !i %in% name.extend.par){
          tem <- delta.fun(data.par[par.id[i], 'link'], tem, o_value[['D']])
          names(tem) <- 'D'
          out_se[[i]] <- c(out_se[[i]], tem)
        }
      } else {
        se_tmb[[i]] <- NA
        names(se_tmb[[i]]) <- name_output_tmb[[i]]
      }
    } else if(i == 'esa') {
      out_se[[i]] <- o_sd[[i]]
      names(out_se[[i]]) <- name_output[[i]]
    }
    
  }
  
  out$se <- do.call('c', out_se)
  names(out$se) <- do.call('c', lapply(out_se, function(x) names(x)))
  out[['output.tmb']][['se_link']] <- se_tmb
  #####################################################################################################
  ## Insert the necessary TMB output.
  tem <- vector('list', 2)
  names(tem) <- c('DX_non_mask', 'DX_mask')
  tem[['DX_non_mask']] <- DX.full
  tem[['DX_mask']] <- DX.mask
  
  out[['output.tmb']][['DX_standardized']] <- tem
  
  tem <- vector('list', 2)
  names(tem) <- c('value', 'std')
  tem[['value']] <- val(o_value$esa)
  tem[['std']] <- o_sd$esa
  
  out[['esa']] <- tem
  
  out[['output.tmb']][['param.og']] <- param.og
  out[['output.tmb']][['param.extend']] <- name.extend.par
  out[['output.tmb']][['param.fix']] <- name.fixed.par
  out[['output.tmb']][['param.info.table']] <- data.par
  out[['output.tmb']][['data.traps']] <- data.traps
  out[['output.tmb']][['dims']] <- dims
  out[['output.tmb']][['area_unit']] <- A
  out[['output.tmb']][['survey.length']] <- survey.length
  out[['output.tmb']][['avg_cue_rates']] <- mean.cue.rates
  out[['output.tmb']][['data.full']] <- data.full
  out[['output.tmb']][['data.mask']] <- data.mask
  out[['output.tmb']][['data.dists.thetas']] <- data.dists.thetas
  out[['output.tmb']][['detfn']] <- detfn
  out[['output.tmb']][['ss.link']] <- ss.link
  out[['output.tmb']][['cutoff']] <- cutoff
  out[['output.tmb']][['sound.speed']] <- sound.speed
  out[['output.tmb']][['gam_output']] <- gam.output
  out[['output.tmb']][['tmb_output_og']] <- tmb.output.og
  out[['output.tmb']][['mean_std_for_standardizing']] <- lst.mean.std
  
  ######################################################################################################
  ## The 5th component: "loglik".
  out$loglik <- -1 * opt$objective
  
  out$opt_output <- opt
  
  #####################################################################################################
  ## The 6th component: "maxgrad".
  out$maxgrad <- max(abs(o$gradient.fixed))
  
  #####################################################################################################
  ## The 7th and 8th components: "cor" and "vcov"; after discussion, it is fine to
  ## use the output from TMB directly.
  tem <- name_output_tmb
  
  for(i in names(tem)){
    if(i %in% name.fixed.par) tem[i] <- list(NULL)
  }
  tem <- do.call('c', tem)
  tem <- c(tem, paste('esa', seq(dims$n.sessions), sep = "."))
  ## Remove fixed parameters from the covariance matrix.
  if(length(name.fixed.par) > 0){
    which_fixed <- which(names(o$value) %in% name.fixed.par)
    out$vcov <- o$cov[-which_fixed, -which_fixed]
  } else {
    out$vcov <- o$cov
  }
  
  dimnames(out$vcov) <- list(tem, tem)
  
  ## Convert the covariance matrix to a correlation matrix.
  tem <- sqrt(diag(out$vcov))
  out$cor <- out$vcov / outer(tem, tem)
  
  
  ######################################################################################################
  ## The 9th component: "npar".
  
  out$npar <- sum(sapply(name_output, length))
  
  ######################################################################################################
  ## The 10th to 14th components: not sure what they are, so skip them for now.
  
  
  
  
  ######################################################################################################
  ## The 15th component: 'eratio';
  ## obtain all estimated values excluding 'esa'.
  tem <- o$value[which(!names(o$value) %in% c('esa', name.fixed.par))]
  out$eratio <- max(abs(o$gradient.fixed) / abs(tem))
  
  
  ######################################################################################################
  ## The 16th component: "D.mask".
  out$D.mask <- vector('list', dims$n.sessions)
  
  for(s in 1:dims$n.sessions){
    tem <- paste0("D_mask[", s, "]")
    out$D.mask[[s]] <- out$coeflist[[tem]]
    tem <- stringr::str_pad(as.character(1:dims$n.masks[s]), nchar(as.character(dims$n.masks[s])),
                           pad = '0', side = 'left')
    names(out$D.mask[[s]]) <- paste(paste0('D_mask[', s, ']'), tem, sep = ".")
  }
  
  
  
  ########################################################################################################
  ## The 17th component: "mm.ihd".
  
  out$mm.ihd <- vector('list', dims$n.sessions)
  
  for(s in 1:dims$n.sessions){
    ## First get the first row of this session in data.trap.
    tem <- min(which(data.trap$session == s))
    ## Then get the data from DX.full, the design matrix for data.trap.
    tem <- DX.full$D[tem,]
    ## Finally, generate a matrix with n.masks rows that contains the information
    ## on "D" at the session level.
    tem <- matrix(rep(tem, dims$n.masks[s]), nrow = dims$n.masks[s], byrow = T)
    
    if(data.par[par.id['D'], 'n_col_mask'] > 0){
      out$mm.ihd[[s]] <- cbind(tem, DX.mask$D[data.mask$session == s,])
    } else {
      out$mm.ihd[[s]] <- tem
    }
    colnames(out$mm.ihd[[s]]) <- gsub("D\\.", "", name_output$D)
    rownames(out$mm.ihd[[s]]) <- NULL
  }
  
  
  ######################################################################################################
  ## The 18th component: "args".
  
  out$args <- arg.input
  
  ######################################################################################################
  ## The 19th component: 'n.sessions'.
  out$n.sessions <- dims$n.sessions
  #####################################################################################################
  ## The 20th component: 'fit.types'.
  tem <- logical(5)
  names(tem) <- c('bearing', 'dist', 'ss', 'toa', 'mrds')
  
  tem['bearing'] <- ('bearing' %in% bucket.info)
  tem['dist'] <- ('dist' %in% bucket.info)
  tem['ss'] <- ('ss' %in% bucket.info)
  tem['toa'] <- ('toa' %in% bucket.info)
  tem['mrds'] <- ('mrds' %in% bucket.info)
  
  out$fit.types <- tem
  
  #####################################################################################################
  ## The 21st component: 'infotypes'.
  
  out$infotypes <- bucket.info[which(bucket.info %in% c('bearing', 'dist', 'ss', 'toa', 'mrds'))]
  
  
  #####################################################################################################
  ## The 22nd component: 'detpars'.
  out$detpars <- param.og[which(!param.og %in% c('kappa', 'alpha', 'sigma.toa', 'D'))]
  
  
  ######################################################################################################
  ## The 23rd component: 'suppars'.
  out$suppars <- param.og[which(param.og %in% c('kappa', 'alpha', 'sigma.toa'))]
  if(length(out$suppars) == 0) out$suppars <- NULL
  
  ######################################################################################################
  ## The 24th component: 'D.betapars'.
  out$D.betapars <- name_output$D
  
  ######################################################################################################
  ## The 25th component: 'phases'.
  all.pars <- c(out$detpars, out$suppars, out$D.betapars)
  
  out$phases <- ifelse(all.pars %in% names(out$args$fix), -1, 1)
  names(out$phases) <- all.pars
  out$phases <- as.list(out$phases)
  
  #####################################################################################################
  ## The 26th component: 'par.links';
  ## the 27th component: 'par.unlinks'.
  
  out$par.links <- vector('list', length(all.pars))
  names(out$par.links) <- all.pars
  
  out$par.unlinks <- vector('list', length(all.pars))
  names(out$par.unlinks) <- all.pars
  
  for(i in c(out$detpars, out$suppars)){
    if(!i %in% name.extend.par){
      if(data.par[par.id[i], 'link'] == 'log'){
        out$par.links[[i]] <- cus.log
        out$par.unlinks[[i]] <- cus_log_unlink
      } else if(data.par[par.id[i], 'link'] == 'logit'){
        out$par.links[[i]] <- cus.logit
        out$par.unlinks[[i]] <- cus.logit.unlink
      } else {
        out$par.links[[i]] <- cus.identity
        out$par.unlinks[[i]] <- cus.identity
      }
    } else {
      out$par.link[[i]] <- cus.identity
      out$par.unlinks[[i]] <- cus.identity
    }
  }
  
  for(i in out$D.betapars){
    out$par.links[[i]] <- cus.identity
    out$par.unlinks[[i]] <- cus.identity
  }
  
  
  ######################################################################################################
  ## The 28th component: 'fit.ihd'.
  ## Only if the model is 'D' extended and the extension is a "mask" extension only
  ## (in other words, the model is equivalent to the incumbent "inhomo" model)
  ## is this 'ihd' set to TRUE.
  
  if('D' %in% name.extend.par & data.par[par.id['D'], 'n_col_full'] == 1){
    out$fit.ihd <- TRUE
  } else {
    out$fit.ihd <- FALSE
  }
  
  
  #####################################################################################################
  ## The 29th component: 're.detfn'. Not sure what this is;
  ## in all of the tested sample data it is FALSE, so set it to FALSE for now.
  
  out$re.detfn <- FALSE
  
  
  #####################################################################################################
  ## The 30th component: 'fit.freqs';
  ## this is not included in the TMB template yet.
  
  out$fit.freqs <- !is.null(cue.rates)
  
  ######################################################################################################
  ## The 31st component: 'first.calls';
  ## based on the instructions, set it to FALSE.
  
  out$first.calls <- FALSE
  
  
  
  #######################################################################################################
  ## The 32nd component: 'model.formula'.
  if('D' %in% name.extend.par){
    out$model.formula <- stats::as.formula(paste(c('gam.resp', as.character(par.extend$model$D)), collapse = ''))
  } else {
    out$model.formula <- stats::as.formula('gam.resp~1')
  }
  
  ########################################################################################################
  ## The 33rd component: 'fgam'.
  out$fgam <- fgam
  
  ########################################################################################################
  ## The 34th component: 'all.covariates'.
  if(is(out$args$mask, 'list')){
    out$all.covariates <- do.call('rbind', out$args$mask)
  } else {
    out$all.covariates <- out$args$mask
  }
  
  if(is.scale){
    for(i in 1:ncol(out$all.covariates)){
      tem <- out$all.covariates[, i]
      out$all.covariates[, i] <- (tem - mean(tem)) / sd(tem)
    }
  }
  
  if(ncol(fgam$X) > 1){
    out$all.covariates <- cbind(out$all.covariates, fgam$X[,-1])
  }
  out$all.covariates <- as.data.frame(out$all.covariates)
  
  ###########################################################################################################
  ## The 35th/36th components: conditional likelihood and partial derivatives.
  out$CL <- CL
  out$esa_partial_derivatives <- esa.partial.derivatives
  
  ############################################################################################################
  return(out)
}

#' Buffered trace printer for TMB objectives (paired fn/grad)
#'
#' Returns a small tap object with callbacks that you attach to the objective
#' and gradient produced by [TMB::MakeADFun()]. While optimizing (e.g. via
#' [stats::nlminb()]), it prints a right-aligned, fixed-width row **only when**
#' the objective value and gradient were both evaluated at the **same**
#' parameter vector. This guarantees the `mgc` column (maximum absolute gradient
#' component) corresponds to the printed `fval`.
#'
#' Because optimizers often evaluate the objective and gradient at different
#' points and in different orders (e.g. line search), the tap uses an internal
#' buffer to: (a) remember the most recent `fn` and `gr` evaluations, and
#' (b) print a row only when those evaluations match in parameters (within
#' `mgc.tol`).
#'
#' Columns are `iter` (count of printed rows), `fval`, optional `mgc`, optional
#' `step` (Euclidean or max coordinate change since the previous *printed* row),
#' followed by the parameter columns in `trace.cols`.
#'
#' @param trace.cols Character vector of parameter names to print as columns
#'   (typically `names(obj$par)` or a subset).
#' @param only.improvements Logical; if `TRUE` (default) print only when `fval`
#'   strictly improves over the best so far (an “outer-like” trace). Set `FALSE`
#'   to print every paired evaluation.
#' @param step One of `"euclid"`, `"max"`, `"none"`. Controls the optional
#'   step-size column between successive *printed* parameter vectors.
#' @param show.mgc Logical; if `TRUE`, include a `mgc` column equal to
#'   `max(abs(grad))` at the matched point. Set `FALSE` if you do not have a
#'   gradient (e.g. `gr.skip = TRUE`).
#' @param mgc.tol Numeric tolerance used to decide whether `fn` and `gr`
#'   parameter vectors “match” (default `1e-10`).
#'
#' @return A list with callbacks:
#' \itemize{
#'   \item `on_fn(par_vec, fval, par_names, G_matrix_list)`: 
#'   record an objective evaluation.
#'   \item `on_gr(par_vec, grad_vec, par_names, G_matrix_list)`: 
#'   record a gradient evaluation.
#' }
#' Each callback returns `invisible(NULL)`; printing happens as a side effect.
#'
#' @section Usage:
#' \preformatted{
#' tap <- make.buffer.printer(trace.cols = names(obj$par),
#'                            step = "euclid",
#'                            show.mgc = TRUE)
#' obj$fn <- tap.fn(obj$fn, tap, par.names = names(obj$par))
#' obj$gr <- tap.gr(obj$gr, tap, par.names = names(obj$par))
#'
#' opt <- stats::nlminb(obj$par, obj$fn, obj$gr, 
#'                                control = list(trace = 0))
#' }
#'
#' @notes
#' - `iter` counts **printed rows** (accepted, paired evaluations), not
#'   `nlminb`'s internal major-iteration counter. Official counts are in
#'   `opt$iterations` / `opt$evaluations`.
#' - `mgc` and `step` are on the **TMB parameterization** scale.
#'
#' @seealso [TMB::MakeADFun()], [stats::nlminb()]
#' @keywords internal
make.buffer.printer <- function(trace.cols,
                                only.improvements = TRUE,
                                step = c("euclid","max","none"),
                                show.mgc = TRUE,
                                mgc.tol = 1e-10) {
  ## Method to be used when calculating step size.
  step <- match.arg(step)
  
  ## Initialize the printer state.
  st <- new.env(parent = emptyenv())
  st$printed_header <- FALSE
  st$best_f <- Inf
  st$iter_print <- 0L
  st$prev_print_par <- NULL
  st$busy <- FALSE
  
  ## Fn and gr buffers;
  ## required in order to match the corresponding gradient calculation to the parameter
  ## vector.
  st$buf_fn_par <- NULL
  st$buf_fn_f <- NULL
  st$buf_gr_par <- NULL
  st$buf_gr_mgc <- NA_real_
  
  ## Determine columns to be printed.
  cols <- c("iter","fval")
  if (isTRUE(show.mgc)) {
    cols <- append(cols, "mgc", after = 2L)
  }
  if (step != "none") {
    cols <- append(cols, "step", after = match("fval", cols))
  }
  ## Add columns for each parameter.
  cols <- c(cols, trace.cols)
  
  ## Base formatting:
  ## iter - %6d - integer in a field width of 6
  ## fval - %12.6f - fixed-point number, width 12 with 6 decimals
  ## rest - %10.4g - general format, 4 significant digits in a field width of 10.
  fmts <- c(iter = "%6d", fval = "%12.6f", mgc = "%10.4g", step = "%10.4g",
            setNames(rep("%10.4g", length(trace.cols)), trace.cols))
  
  ## Pulls the width out of a format of similar form to "%6d" above.
  minw <- function(fmt) {
    w <- suppressWarnings(as.integer(sub(".*%([0-9]+).*", "\\1", fmt)))
    ## If there are no digits immediately after "%", w will be NA.
    if (is.na(w)) 0L else w
  }
  
  ## Compute the final column widths;
  ## this ensures that the printed columns have the same width as the column headers.
  widths <- mapply(function(h, f) max(nchar(h), minw(f)), cols, fmts[cols])
  
  ## Rebuild fmts so their width matches the computed width (keep precision/type).
  mkfmt <- function(fmt, width) {
    ## Regex that reconstructs width formatting like:
    ## % [flags] [width] [.prec] [type].
    sub("^%([-+ 0#]*)([0-9]*)(\\.?[0-9]*)([a-zA-Z])$",
        paste0("%\\1", width, "\\3\\4"),
        fmt)
  }
  ## Update with the new formatting.
  fmts_aligned <- mapply(mkfmt, fmts[cols], widths, SIMPLIFY = TRUE, USE.NAMES = TRUE)
  
  ## The header will be printed when the `printed_header` flag is FALSE.
  header.print <- function() {
    pieces <- mapply(function(h, w) format(h, width=w, justify="right"),
                     cols, widths, SIMPLIFY=TRUE)
    cat(paste(pieces, collapse=" "), "\n")
    ## `flush.console()` is here to make sure printing still works if the user
    ## is working in the console. There is potential for odd print formatting if it is not
    ## included.
    utils::flush.console()
  }
  
  ## Euclid corresponds to the Euclidean norm, and max to the maximum absolute difference
  ## over the parameters being optimised.
  step_size <- switch(step,
                      euclid = function(a,b) if (is.null(b)) NA_real_ else sqrt(sum((a-b)^2)),
                      max = function(a,b) if (is.null(b)) NA_real_ else max(abs(a-b)),
                      none = function(a,b) NA_real_
  )
  
  ## The actual print function.
  maybe.print <- function(par.names) {
    ## If there is no fn buffer currently, do nothing.
    if (is.null(st$buf_fn_par)) return(invisible())
    
    ## If showing mgc, require a matching gr buffer.
    if (isTRUE(show.mgc)) {
      if (is.null(st$buf_gr_par)) return(invisible())
      ## Check if the fn buffer matches the gr buffer by comparing the
      ## current parameter estimate values, and making sure the difference
      ## between them is negligible.
      same <- length(st$buf_fn_par) == length(st$buf_gr_par) &&
        max(abs(st$buf_fn_par - st$buf_gr_par)) <= mgc.tol
      ## If they are different, again do nothing.
      if (!same) return(invisible())
    }
    
    ## Retrieve the current function value (negative log-likelihood).
    fval <- st$buf_fn_f
    if (only.improvements) {
      ## If the NLL is infinite (edge case), or if it is less than the best
      ## NLL so far, then ignore this step and reset the buffers.
      if (!is.finite(fval) || fval >= st$best_f - 1e-12) {
        st$buf_fn_par <- NULL 
        st$buf_fn_f <- NULL
        ## Move on to the next step.
        return(invisible())
      }
    }
    st$best_f <- min(st$best_f, fval)
    
    ## Print the header if it hasn't been printed yet.
    if (!st$printed_header) { 
      header.print()
      ## Update the flag so the header is not printed again.
      st$printed_header <- TRUE 
    }
    st$iter_print <- st$iter_print + 1L
    
    ## Build the row to be printed.
    row <- setNames(rep(NA_real_, length(cols)), cols)
    row["iter"] <- st$iter_print
    row["fval"] <- fval
    if ("mgc" %in% cols) {
      row["mgc"] <- st$buf_gr_mgc
    }
    if ("step" %in% cols) {
      row["step"] <- step_size(st$buf_fn_par, st$prev_print_par)
    }
    st$prev_print_par <- st$buf_fn_par
    
    keep <- intersect(trace.cols, par.names)
    if (length(keep)) {
      row[keep] <- st$buf_fn_par[match(keep, par.names)]
    }
    ## Build the final row string.
    pieces <- mapply(function(v, fmt) sprintf(fmt, v),
                     row, fmts_aligned, SIMPLIFY = TRUE)
    ## Print the row.
    cat(paste(pieces, collapse=" "), "\n")
    utils::flush.console()
    
    ## Clear buffers.
    st$buf_fn_par <- NULL
    st$buf_fn_f <- NULL
    st$buf_gr_par <- NULL
    st$buf_gr_mgc <- NA_real_
    invisible()
  }
  
  on.fn.evaluation <- function(par.vec, fval, par.names) {
    ## If there is already a connection to the buffer environment,
    ## then make sure recursive looping is avoided.
    if (st$busy) return(invisible(NULL))
    st$busy <- TRUE
    on.exit({ st$busy <- FALSE }, add = TRUE)
    
    ## Fill the function-side buffer; the gradient side may arrive earlier or later.
    st$buf_fn_par <- par.vec
    st$buf_fn_f <- fval
    
    ## Attempt to print; this only prints when the required buffers are present
    ## and (if `show.mgc=TRUE`) par vectors match within tolerance.
    maybe.print(par.names)
    invisible(NULL)
  }
  
  on.gr.evaluation <- function(par.vec, grad.vec, par.names) {
    ## See above `on.fn.evaluation()`.
    if (st$busy) return(invisible(NULL))
    st$busy <- TRUE
    on.exit({ st$busy <- FALSE }, add = TRUE)
    st$buf_gr_par <- par.vec
    st$buf_gr_mgc <- max(abs(grad.vec))
    maybe.print(par.names)
    invisible(NULL)
  }
  
  list(on_fn_evaluation = on.fn.evaluation, on_gr_evaluation = on.gr.evaluation)
  
}

#' Wrap an objective to emit buffered trace rows on each evaluation
#'
#' Creates an idempotent wrapper around an objective function `fn(par, ...)`
#' that forwards the call unchanged, then notifies a tap object produced by
#' [make.buffer.printer()] via `tap$on_fn_evaluation(par_vec, fval, par.names)`.
#' Use this with [stats::nlminb()] (or similar optimizers) to stream live,
#' correctly paired trace rows when combined with [tap.gr()].
#'
#' The wrapper is **idempotent**: if it detects the function was already wrapped
#' (attribute `.__fn_tapped`), it returns `fn` unchanged. Arguments are
#' `force()`d so the wrapper closes over the current objects.
#'
#' @param fn Objective function with signature `function(par, ...) -> scalar`.
#' @param tap Tap/list returned by [make.buffer.printer()], expected to provide
#'   a method `on_fn_evaluation(par_vec, fval, par.names)`.
#' @param par.names Character vector of names corresponding to `par`. These are
#'   used by the tap to map positions in `par` onto the `trace_cols` it prints.
#' @param G.matrix.list List of back-transformation matrices used to back 
#'   transform parameter values in the case they have been scaled and centered
#'   for model fitting.
#'
#' @return A function with the same signature and return value as `fn`, which
#'   also triggers the tap side-effect after each evaluation.
#'
#' @seealso [tap.gr()], [make.buffer.printer()], [stats::nlminb()]
#' @keywords internal
tap.fn <- function(fn, tap, par.names, G.matrix.list) {
  ## Make sure we evaluate the arguments as they are currently.
  force(fn)
  force(tap)
  force(par.names)
  force(G.matrix.list)
  
  ## Make sure that the function is not already wrapped.
  if (isTRUE(attr(fn, ".__fn_tapped"))) return(fn)
  
  ## Wrap the function.
  wrapped <- function(p, ...) {
    f <- fn(p, ...)
    
    u_names <- unique(names(p))
    for(i in u_names) {
      index_u_name <- which(names(p) == i)
      p[index_u_name] <- G.matrix.list[[i]] %*% p[index_u_name]
    }
    tap$on_fn_evaluation(as.numeric(p), f, par.names)
    f
  }
  
  attr(wrapped, ".__fn_tapped") <- TRUE
  wrapped
}

#' Wrap a gradient to emit buffered trace rows (mgc) on each evaluation
#'
#' Creates an idempotent wrapper around a gradient function `gr(par, ...)`
#' that forwards the call unchanged, then notifies a tap object produced by
#' [make.buffer.printer()] via `tap$on_gr_evaluation(par_vec, grad_vec, par.names)`.
#' When used together with [tap.fn()], the tap can pair `fn` and `gr`
#' evaluations at the same parameter vector and print a trace row that includes
#' `mgc = max(abs(grad))`.
#'
#' If `gr` is `NULL`, this helper returns `NULL` (no wrapping).
#' The wrapper is **idempotent**: if it detects the function was already wrapped
#' (attribute `.__gr_tapped`), it returns `gr` unchanged. Arguments are
#' `force()`d so the wrapper closes over the current objects.
#'
#' @param gr Gradient function with signature `function(par, ...) -> numeric`
#'   (same length/order as `par`), or `NULL` if no gradient is supplied.
#' @param tap Tap/list returned by [make.buffer.printer()], expected to provide
#'   a method `on_gr_evaluation(par_vec, grad_vec, par.names)`.
#' @param par.names Character vector of names corresponding to `par`. These are
#'   used by the tap to map gradient components to parameter columns.
#' @param G.matrix.list List of back-transformation matrices used to back 
#'   transform parameter values in the case they have been scaled and centered
#'   for model fitting.
#'
#' @return A function with the same signature and return value as `gr`,
#'   which also triggers the tap side-effect after each evaluation; or `NULL`
#'   if `gr` was `NULL`.
#'
#' @seealso [tap.fn()], [make.buffer.printer()], [stats::nlminb()]
#' @keywords internal
tap.gr <- function(gr, tap, par.names, G.matrix.list) {
  ## Make sure we evaluate the arguments as they are currently.
  if (is.null(gr)) return(NULL)
  force(gr)
  force(tap)
  force(par.names)
  force(G.matrix.list)
  
  ## Make sure that the function is not already wrapped.
  if (isTRUE(attr(gr, ".__gr_tapped"))) return(gr)
  
  ## Wrap the function.
  wrapped <- function(p, ...) {
    g <- gr(p, ...)
    
    u_names <- unique(names(p))
    for(i in u_names) {
      index_u_name <- which(names(p) == i)
      p[index_u_name] <- G.matrix.list[[i]] %*% p[index_u_name]
    }
    
    tap$on_gr_evaluation(as.numeric(p), g, par.names)
    g
  }
  
  attr(wrapped, ".__gr_tapped") <- TRUE
  wrapped
}






