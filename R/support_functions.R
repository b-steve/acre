fulllist.par.generator <- function(){
  ## The order matters and cannot be changed.
  return(c('g0', 'sigma', 'lambda0', 'z', 'shape.1',
            'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
            'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa',
            "sigma.b0.ss", 'D', 'mu'))
}

natural.number.check <- function(session, key){
  output <- TRUE
  for(s in unique(session)){
    tem_key <- key[which(session == s)]
    seq_key <- seq(length(unique(tem_key)))
    if(any(!tem_key %in% seq_key)) output <- FALSE
  }

  return(output)
}

convert.natural.number <- function(dat, is.animalID, which.convert = 'both'){

  dat.na <- subset(dat, is.na(dat$ID))
  dat <- subset(dat, !is.na(dat$ID))

  if(which.convert == 'ID'){
    if(is.animalID){
      keys <- paste(dat$session, dat$animal_ID, sep = '_')
    } else {
      keys <- dat$session
    }
    uni_keys <- unique(keys)
    output <- vector('list', length(uni_keys))
    for(i in 1:length(uni_keys)){
      output[[i]] <- dat[keys == uni_keys[i],,drop = FALSE]
      output[[i]]$ID <- as.numeric(as.factor(output[[i]]$ID))
    }
    output <- do.call('rbind', output)
  } else if(which.convert == 'animal_ID'){
    stopifnot(is.animalID)
    keys <- dat$session
    uni_keys <- unique(keys)
    output <- vector('list', length(uni_keys))
    for(i in 1:length(uni_keys)){
      output[[i]] <- dat[keys == uni_keys[i],,drop = FALSE]
      output[[i]]$animal_ID <- as.numeric(as.factor(output[[i]]$animal_ID))
    }
    output <- do.call('rbind', output)
  } else if(which.convert == 'both'){
    if(is.animalID){
      dat <- convert.natural.number(dat, is.animalID, 'animal_ID')
    }
    dat <- convert.natural.number(dat, is.animalID, 'ID')
    output <- dat
  } else {
    stop('Invalid input.')
  }

  ## Rows from sessions without detections have NA 'ID' (and 'animal_ID'), but the column may still be
  ## character; if so, rbind() would turn the converted IDs back to character, which are then sorted
  ## lexicographically ("1", "10", "2", ...) and no longer match the numeric ID used to index data_full in TMB.
  if(which.convert %in% c('ID', 'both')) dat.na$ID <- as.numeric(dat.na$ID)
  if(is.animalID & which.convert %in% c('animal_ID', 'both')) dat.na$animal_ID <- as.numeric(dat.na$animal_ID)

  output <- rbind(output, dat.na)
  return(output)

}


agg.sort <- function(dat, obj, lst, f){
  if(length(obj) == 1){
    x <- dat[[obj]]
  } else {
    x <- dat[[obj[1]]]
    for(i in 2:length(obj)){
      x <- paste(x, dat[[obj[i]]], sep = '---')
    }
  }

  y <- vector('list', length(lst))
  names(y) <- lst
  for(i in lst) y[[i]] <- dat[[i]]
  output <- aggregate(x, y, function(x) f(x))

  sort_by <- vector('list', length(lst))
  for(i in 1:length(lst)) sort_by[[i]] <- output[[lst[i]]]

  o <- do.call('order', sort_by)
  output <- output[o,]
  return(output)
}


df.to.list <- function(df, n.sessions){
  if(is(df, 'list')){
    stopifnot(length(df) == n.sessions)
    return(df)
  } else {
    stopifnot(any(is(df, 'data.frame'), is(df, 'matrix')))
    output <- vector('list', n.sessions)
    for(s in 1:n.sessions) output[[s]] <- df
    return(output)
  }
}


extend.dat.check <- function(dat, check.var, ori.dat, n.sessions, n.var, identical.flag){
  stopifnot(is(dat, 'data.frame'))
  if(!(check.var %in% colnames(dat))){
    stopifnot(nrow(dat) == n.var[1])
    dat[[check.var]] <- seq(n.var[1])
  }

  if('session' %in% colnames(dat)){
    stopifnot(length(colnames(dat)) > 2)
  } else {
    ## If all sessions share the same traps/mask, the user can omit 'session' from the traps/mask covariates,
    ## or there is only 1 session.
    stopifnot(length(colnames(dat)) > 1)
    if(n.sessions > 1) stopifnot(identical.flag)
    tem <- vector('list', n.sessions)
    for(s in 1:n.sessions){
      tem[[s]] <- dat
      tem[[s]][['session']] <- s
    }
    dat <- do.call('rbind', tem)
  }

  if(any(duplicated(dat[, c('session', check.var)]))) stop(paste0('duplicated ', check.var, ' input.'))

  stopifnot(all(paste(dat$session, dat[[check.var]], sep = '-') %in%
                  unique(paste(ori.dat$session, ori.dat[[check.var]], sep = '-'))))

  return(dat)
}



covariates.mask.check <- function(dat, n.sessions, n.masks, identical.flag){
  stopifnot(is(dat, 'list') | is(dat, 'data.frame'))

  if(is(dat, 'list')){
    if(length(dat) > 1){
      ## If there is more than one component in this list, each component represents one session.
      stopifnot(length(dat) == n.sessions)
      stopifnot(all(sapply(dat, nrow) == n.masks))
    } else {
      ## If there is only one component in this list, then n.sessions > 1 is a special scenario.
      stopifnot(nrow(dat[[1]]) == n.masks[1])

      if(n.sessions > 1){
        stopifnot(identical.flag)
        tem <- dat[[1]]
        dat <- vector('list', n.sessions)
        for(s in 1:n.sessions) dat[[s]] <- tem
      }

    }

  } else {
    stopifnot(nrow(dat) == n.masks[1])
    if(n.sessions > 1){
      stopifnot(identical.flag)
    }

    tem <- dat
    dat <- vector('list', n.sessions)
    for(s in 1:n.sessions) dat[[s]] <- tem
  }

  for(s in 1:n.sessions){
    dat[[s]][['session']] <- s
    dat[[s]][['mask']] <- 1:n.masks[s]
  }
  dat <- do.call('rbind', dat)

  return(dat)
}


default.link <- function(i){
  if(i == "mu") return('log')
  if(i == "D") return('log')
  if(i == "g0") return('logit')
  if(i == "sigma") return('log')
  if(i == "lambda0") return('log')
  if(i == "z") return('log')
  if(i == "shape.1") return('log')
  if(i == "shape.2") return('identity')
  if(i == "shape") return('identity')
  if(i == "scale") return('log')
  if(i == "b0.ss") return('log')
  if(i == "b1.ss") return('log')
  if(i == "b2.ss") return('log')
  if(i == "sigma.ss") return('log')
  if(i == "kappa") return('log')
  if(i == "alpha") return('log')
  if(i == "sigma.toa") return('log')
  if(i == "sigma.b0.ss") return('log')
}



detfn.params <- function(detfn){
  if(detfn == 'hn'){
    param.og <- c('g0', 'sigma')
  }
  if(detfn == 'hhn'){
    param.og <- c('lambda0', 'sigma')
  }
  if(detfn == 'hr'){
    param.og <- c('g0', 'sigma', 'z')
  }
  if(detfn == 'th'){
    param.og <- c('shape', 'scale')
  }
  if(detfn == 'lth'){
    param.og <- c('shape.1', 'shape.2', 'scale')
  }
  if(detfn == 'ss'){
    param.og <- c('b0.ss', 'b1.ss', 'sigma.ss')
  }

  return(param.og)
}

## See section 1.4.2 in:
## https://research-repository.st-andrews.ac.uk/bitstream/handle/10023/18233/BenStevensonPhDThesis.pdf.
p.dot.defaultD <- function(points = NULL, traps = NULL, detfn = NULL,
                          ss.link = NULL, pars = NULL, A, n.quadpoints = 8,
                          esa = T) {
  dists <- distances(traps, points)
  probs <- det.prob(detfn, pars, dists, ss.link)
  if(detfn == 'ss'){
    sigma.ss <- pars$sigma.ss
    cutoff <- pars$cutoff
    probs <- 1 - pnorm(cutoff, mean = probs, sd = sigma.ss)
  }
  ## The term prod(1 - x) is P(didn't detect on trap 1 & didn't detect on trap 2 & ...), and
  ## 1 - prod(1 - x) is P(detect on at least 1 trap).
  out <- 1 - apply(1 - probs, 2, prod)

  ## If esa == TRUE: returns the ESA.
  ## If esa == FALSE: returns P(detect on at least 1 trap) for each mask point.
  if (esa) {
    out <- A*sum(out)
    return(out)
  } else {
    return(out)
  }
}

#' Calulating probabilities of detection at mask points
#'
#' Calculates the probability of a call produced at a mask point being
#' detected by one or more detectors.
#' @inheritParams boot.acre
#'
#' @return A list with a component for each session, containing a
#'   vector of detection probabilities.
#'
#' @export
p.dot <- function(fit){
  n.sessions <- fit$n.sessions
  out <- vector(mode = "list", length = n.sessions)
  for (i in 1:n.sessions){
    mask <- as.matrix(get.mask(fit)[[i]])
    traps <- get.trap(fit)[[i]]
    det.pars <- as.list(coef(fit, type = "fitted"))
    ext.par.matrices <- get.par.extend.matrix(fit, mask, traps, session = i)
    det.pars[names(ext.par.matrices)] <- ext.par.matrices
    out[[i]] <- p.dot.defaultD(points = mask, traps = traps,
                               detfn = get.detfn(fit),
                               ss.link = get.ss.link(fit),
                               pars = det.pars,
                               A = attr(mask, "area"),
                               esa = FALSE)
  }
  out
}

separate.formula <- function(foo, var.m){
  foo_vars <- all.vars(foo[[3]])
  response_var <- all.vars(foo[[2]])
  index <- which(foo_vars %in% var.m)
  if(length(index) == 0){
    full_vars <- foo_vars
    mask_vars <- character(0)
  } else {
    full_vars <- foo_vars[-index]
    mask_vars <- foo_vars[index]
  }

  foo_terms <- attr(stats::terms(foo), 'factors')[-1,,drop = FALSE]

  row_names <- as.list(rownames(foo_terms))
  col_names <- as.list(colnames(foo_terms))
  ## Use column names to check whether there is an interaction term
  ## between variables from mask level and other levels.
  for(i in 1:length(col_names)){
    tem.name <- col_names[[i]]
    tem.foo <- stats::as.formula(paste0('gam.resp~',tem.name))
    tem.foo_vars <- all.vars(tem.foo[[3]])
    if(!(sum(tem.foo_vars %in% var.m) %in% c(0, length(tem.foo_vars)))){
      stop("interaction between variables from mask level and other level is not supported.")
    }
  }
  ## Use row names to separate this foo_terms matrix into two matrices,
  ## one for the mask level and the other for the other levels.
  for(i in 1:length(row_names)){
    tem.name <- row_names[[i]]
    tem.foo <- stats::as.formula(paste0('gam.resp~',tem.name))
    tem.foo_vars <- all.vars(tem.foo[[3]])
    row_names[[i]] <- tem.foo_vars
  }

  len <- lapply(row_names, length)
  foo_terms <- foo_terms[rep(1:nrow(foo_terms), len),,drop = FALSE]
  row_names <- do.call('c', row_names)
  rownames(foo_terms) <- row_names

  n.var.full <- sum(!(foo_vars %in% var.m))
  n.var.mask <- sum((foo_vars %in% var.m))

  if(n.var.full > 0) {
    tem_full <- vector('list', n.var.full)
    names(tem_full) <- full_vars
  }

  if(n.var.mask > 0) {
    tem_mask <- vector('list', n.var.mask)
    names(tem_mask) <- mask_vars
  }
  for(i in foo_vars){
    if(i %in% full_vars){
      tem_full[[i]] <- foo_terms[which(row_names == i),,drop = FALSE]
      tem_full[[i]] <- apply(tem_full[[i]], 2, function(x) sum(x) > 0)
      tem_full[[i]] <- matrix(tem_full[[i]], ncol = ncol(foo_terms),
                             byrow = T)
      colnames(tem_full[[i]]) <- colnames(foo_terms)
      rownames(tem_full[[i]]) <- i
    } else {
      tem_mask[[i]] <- foo_terms[which(row_names == i),,drop = FALSE]
      tem_mask[[i]] <- apply(tem_mask[[i]], 2, function(x) sum(x) > 0)
      tem_mask[[i]] <- matrix(tem_mask[[i]], ncol = ncol(foo_terms),
                             byrow = T)
      colnames(tem_mask[[i]]) <- colnames(foo_terms)
      rownames(tem_mask[[i]]) <- i
    }
  }

  if(n.var.full > 0){
    tem <- apply(do.call('rbind', tem_full), 2, any)
    tem.terms <- names(tem)[which(tem)]
    foo_full <- stats::as.formula(paste0(response_var, "~", paste(tem.terms, collapse = "+")))
  } else {
    foo_full <- stats::as.formula(paste0(response_var, "~ 1"))
  }

  if(n.var.mask > 0){
    tem <- apply(do.call('rbind', tem_mask), 2, any)
    tem.terms <- names(tem)[which(tem)]
    foo_mask <- stats::as.formula(paste0(response_var, "~", paste(tem.terms, collapse = "+")))
  } else {
    foo_mask <- stats::as.formula(paste0(response_var, "~ 1"))
  }

  return(list(foo_full = foo_full, foo_mask = foo_mask))
}


## If the value is valid, return TRUE.

param.range.validate <- function(param, value){
  if(param %in% c("sigma", "lambda0", "z", "shape.1","scale", "sigma.b0.ss", "kappa", "alpha", "sigma.toa",
                  "b0.ss", "b1.ss", "b2.ss", "sigma.ss", 'D', 'mu')){
    if(value >= 0){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if(param == 'g0'){
    if(value >= 0 & value <= 1){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}


## Default start values; this part needs a careful reading of the original function.

default.sv <- function(param, info){
  data.full <- info[['data.full']]
  data.mask <- info[['data.mask']]
  buffer <- info[['buffer']]
  A <- info[['A']]
  ss.opts <- info[['ss.opts']]
  cutoff <- ss.opts[['cutoff']]
  ss.link <- ss.opts[["ss.link"]]
  detfn <- info[['detfn']]
  dims <- info[["dims"]]
  survey.length <- info[["survey.length"]]
  param.og <- info[['param.og']]
  sv <- info[['sv']]
  animal.model <- info[['animal.model']]


  if(param == 'sigma'){
    easy.out <- mean(buffer) / 4
    is.animal_ID <- "animal_ID" %in% colnames(data.full)
    same.traplocs <- (nrow(data.full[!duplicated(data.full[,c('trap_x', 'trap_y')]), ]) == 1)
    if(same.traplocs){
      return(easy.out)
    } else {
      tem <- data.full[!duplicated(data.full[,c("session", "animal_ID"[is.animal_ID], "ID", "trap")]), ]
      tem <- subset(tem, !is.na(tem$bincapt) & tem$bincapt == 1)
      u.id <- unique(tem[, c('session', "animal_ID"[is.animal_ID], "ID")])
      u.id <- paste(u.id$session, u.id$animal_ID[is.animal_ID], u.id$ID, sep = "_")

      mean.dists <- data.frame(avg_dist = numeric(0), weight = numeric(0))
      j <- 1
      for(i in 1:length(u.id)){
        tem1 <- subset(tem, paste(tem$session, tem$animal_ID[is.animal_ID], tem$ID, sep = "_") == u.id[i],
                      select = c("trap_x", "trap_y"))

        if(nrow(tem1) > 1){
          rc.dists <- distances(as.matrix(tem1), as.matrix(tem1))
          mean.dists[j, 1] <- mean(rc.dists[rc.dists > 0])
          mean.dists[j, 2] <- length(rc.dists[rc.dists > 0])
          j <- j + 1
        }
      }

      if(nrow(mean.dists) > 0){
        out <- sum(mean.dists$avg_dist * mean.dists$weight)/sum(mean.dists$weight)
      } else {
        out <- easy.out
      }
      return(out)
    }
  } else if(param == "D"){
    session_to_use <- which(dims$n.IDs > 0)[1]
    tem.mask <- subset(data.mask, data.mask$session == session_to_use)
    mask <- as.matrix(unique(tem.mask[, c('x', 'y')]))
    tem.trap <- subset(data.full, data.full$session == session_to_use)
    traps <- tem.trap[, c('trap', 'trap_x', 'trap_y')]
    traps <- traps[!duplicated(traps[['trap']]), ]
    traps <- as.matrix(traps[order(traps$trap), c('trap_x', 'trap_y')])
    survey.length <- survey.length[session_to_use]
    A <- A[session_to_use]

    detpar.names <- param.og[which(!(param.og %in% c('kappa', 'alpha', 'sigma.toa', 'D')))]
    pars <- sv[detpar.names]



    if(detfn == 'ss'){
      detfn <- 'ss'
      pars[["sigma.b0.ss"]] <- 0
      if(!any(detpar.names == 'b2.ss')) pars[['b2.ss']] <- 0
    }

    if (!is.null(cutoff)){
      pars$cutoff <- cutoff
    }
    esa <- p.dot.defaultD(points = mask, traps = traps, detfn = detfn,
                          ss.link = ss.link, pars = pars, A = A, n.quadpoints = 8, esa = T)
    if(animal.model){
      return(dims$n.animals[session_to_use]/(esa*survey.length))
    } else {
      return(dims$n.IDs[session_to_use]/(esa*survey.length))
    }

    ## This part of the original function is difficult; better to ask for help.
  } else if(param == "g0"){
    return(0.95)
  } else if(param == "lambda0"){
    return(2)
  } else if(param == "z"){
    return(1)
  } else if(param == "sigma.toa"){
    ## Matrix where columns correspond to traps and rows correspond to call IDs.
    toa_by_call_matrix <- matrix(data.full$toa, ncol = max(data.full$trap), byrow = T)
    ## Take the average of the standard deviations of all calls.
    mean(apply(toa_by_call_matrix, 1, sd))
    return(mean(apply(toa_by_call_matrix, 1, sd)))
  } else if(param == "kappa"){
    return(10)
  } else if(param == "b0.ss"){
    if(ss.link == 'log'){
      return(log(max(data.full[['ss']], na.rm = TRUE)))
    } else {
      return(quantile(data.full[['ss']], 0.90, na.rm = TRUE))
    }
  } else if(param == "sigma.b0.ss"){
    if(ss.link == 'log'){
      return(0.01 * log(max(data.full[['ss']], na.rm = TRUE)))
    } else {
      return(0.01 * max(data.full[['ss']], na.rm = TRUE))
    }
  } else if(param == "b1.ss"){
    quantile.ss <- quantile(data.full[['ss']], 0.90, na.rm = TRUE)
    out <- (quantile.ss - cutoff)/(mean(buffer)/2)
    if(ss.link == "log"){
      max.ss <- max(data.full[["ss"]], na.rm = TRUE)
      out <- (max.ss - cutoff)/(mean(buffer)/2)
      out <- out/max.ss
    }
    return(out)
  } else if(param == "b2.ss"){
    return(0.1)
  } else if(param == "sigma.ss"){
    ss_valid <- data.full[['ss']][data.full[['ss']] >= cutoff & !is.na(data.full[['ss']])]
    return(stats::sd(ss_valid))
  } else if(param == "alpha"){
    return(2)
  } else if(param == "shape"){
    return(2)
  } else if(param == "scale"){
    if(detfn == "th"){
      out <- default.sv(param = "sigma", info = info)
    } else if(detfn == "lth"){
      sigma <- default.sv(param = "sigma", info = info)
      shape.1 <- default.sv(param = "shape.1", info = info)
      out <- (log(shape.1) - log(shape.1 - 0.5)) / sigma
    } else {
      stop("Detection function not recognised.")
    }
    return(out)
  } else if(param == "shape.1"){
    return(0.809017)
  } else if(param == "shape.2"){
    sigma <- default.sv(param = "sigma", info = info)
    shape.1 <- default.sv(param = "shape.1", info = info)
    scale <- default.sv(param = "scale", info = info)
    return(shape.1 + scale * sigma)
  } else if(param == 'mu'){
    n.animals <- dims$n.animals
    n.animal.call <- dims$n.animal.call
    n.sessions <- dims$n.sessions

    i <- 1
    for(s in 1:n.sessions){
      if(n.animals[s] > 0){
        for(a in 1:n.animals[s]){
          n.animal.call[i] <- n.animal.call[i]/survey.length[s]
          i <- i + 1
        }
      }
    }

    return(mean(n.animal.call))
  }
}


## Default bounds;
## the default bounds of "D" are not specified in the original.
default.bounds <- function(param){
  if(param == 'g0') {
    return(c(0, 1))
  } else if(param %in% c('sigma', 'lambda0', 'shape.1', 'scale', 'b0.ss',
                         'b1.ss', 'b2.ss', 'sigma.b0.ss', 'sigma.ss',
                         'z', 'sigma.toa', 'alpha', 'mu')){
    return(c(0, 1e8))
  } else if(param == 'kappa'){
    return(c(0, 700))
  } else if(param %in% c('shape', 'shape.2')){
    return(c(-100, 100))
  } else if(param == "D"){
    return(c(0, 1e20))
  }
}

## Link function.
link.fun <- function(link, value){
  if(!link %in% c("identity", "log", "logit")){
    stop('Not a valid link function is assigned.')
  }

  if(link == "identity") return(cus.identity(value))
  if(link == "log") return(cus.log(value))
  if(link == "logit") return(cus.logit(value))

}

unlink.fun <- function(link, value){
  if(!link %in% c("identity", "log", "logit")){
    stop('Not a valid link function is assigned.')
  }

  if(link == "identity") return(value)
  if(link == "log") return(cus_log_unlink(value))
  if(link == "logit") return(cus.logit.unlink(value))
}

delta.fun <- function(link, sd, est){
  if(link == "identity") return(sd)

  ## For link == 'log', we are calculating std(exp(x)) with x_hat = est and x_sd = sd.
  if(link == "log") return(sd * exp(est))

  ## For link == 'logit', we are calculating std(exp(x)/(exp(x) + 1)).
  if(link == 'logit'){
    fx_grad <- exp(est)/(exp(est) + 1) ^ 2
    return(sd * fx_grad)
  }
}

data.sort <- function(dat,  name){
  is.animal_ID <- "animal_ID" %in% colnames(dat)

  if(name == 'data.full'){
    if(is.animal_ID){
      dat <- dat[order(dat$session, dat$animal_ID,
                      dat$ID, dat$trap), ]
    } else {
      dat <- dat[order(dat$session, dat$ID, dat$trap), ]
    }
  }

  if(name == "data.dists.thetas"){
    dat <- dat[order(dat$session, dat$mask, dat$trap), ]
  }

  if(name == "data.ID_mask"){
    if(is.animal_ID){
      ## Calls are ordered as in data.full, so the TMB template can index both models' data the same way.
      dat <- dat[order(dat$session, dat$animal_ID, dat$ID, dat$mask),]
    } else {
      dat <- dat[order(dat$session, dat$ID, dat$mask), ]
    }
  }

  if(name == "data.mask"){
    dat <- dat[order(dat$session, dat$mask),]
  }

  return(dat)
}

numeric.link <- function(link){
  if(is.null(link)){
    ans <- 0
  } else {
    ## If no match, return 0.
    ans <- numeric(length(link))
    for(i in 1:length(link)){
      if(link[i] == 'identity') ans[i] <- 1
      if(link[i] == 'log') ans[i] <- 2
      if(link[i] == 'logit') ans[i] <- 3
      if(link[i] == 'spherical') ans[i] <- 4
    }
  }
  return(ans)
}

numeric.detfn <- function(detfn){
  if(detfn == 'hn') return(1)
  if(detfn == 'hhn') return(2)
  if(detfn == 'hr') return(3)
  if(detfn == 'th') return(4)
  if(detfn == 'lth') return(5)
  if(detfn == 'ss') return(6)
}

cus.log <- function(x){
  x <- pmax(x, .Machine$double.eps)
  return(log(x))
}

cus.logit <- function(x){
  x <- pmax(x, .Machine$double.eps)
  x <- pmin(x, 1 - .Machine$double.eps)
  return(log(x / (1 - x)))
}

cus.identity <- function(x){
  return(x)
}

cus_log_unlink <- .Primitive('exp')

cus.logit.unlink <- function(x){
  x <- pmax(x, cus.logit(0))
  x <- pmin(x, cus.logit(1))
  return(exp(x)/(exp(x) + 1))
}

## Extract values without names from a vector.

val <- function(vec){
  names(vec) <- NULL
  return(vec)
}

## Calculate bearings.
bearings.by.vec <- function(trap.x, trap.y, mask.x, mask.y){
  n <- length(trap.x)
  if(length(trap.y)!=n | length(mask.x)!=n | length(mask.y)!=n){
    stop("different length coordinates.")
  }

  x_diff <- mask.x - trap.x
  y_diff <- mask.y - trap.y

  output <- atan(x_diff / y_diff)

  output <- ifelse(y_diff < 0, output + pi, output)
  output <- ifelse(y_diff >=0 & x_diff < 0, output + 2 * pi, output)

  return(output)
}


## Duplicate the rows k times if a specific column is k,
## in order to split multiple animals/calls (items) from one row into multiple rows;
## used in the simulation function "sim_det_history".
item.split <- function(dat, item){
  ## For each n_call/n_animal > 1, we need to split it into several identical individual calls.
  tem_gt1 <- subset(dat, dat[[item]] > 1)

  if(nrow(tem_gt1) > 0){
    tem_eq1 <- subset(dat, dat[[item]] == 1)

    u_n_items <- unique(tem_gt1[[item]])

    tem_gt1_list <- vector('list', length(u_n_items))

    for(i in 1:length(u_n_items)){
      ## Extract the part where n_items equals this unique number of calls.
      tem <- subset(tem_gt1, tem_gt1[[item]] == u_n_items[i])
      ## Repeat this part u_n_items[i] times.
      tem_list <- vector('list', u_n_items[i])
      for(j in 1:u_n_items[i]) tem_list[[j]] <- tem
      tem_gt1_list[[i]] <- do.call('rbind', tem_list)
    }

    tem_gt1 <- do.call('rbind', tem_gt1_list)

    output <- rbind(tem_gt1, tem_eq1)
    ## Since all calls become individual calls, set the item to 1.
    output[[item]] <- 1
  } else {
    output <- dat
  }

  return(output)
}



list.2vector.4value <- function(param){
  tem <- sapply(param, length)
  output <- do.call('c', param)
  names(output) <- rep(names(tem), tem)
  return(output)
}


ext.par.in.new.df <- function(par.ext, new.covariates, extend.par.covariates){
  output <- NULL
  if(!is.null(par.ext)){
    for(i in par.ext){
      if(all(extend.par.covariates[[i]] %in% colnames(new.covariates))){
        output <- c(output, i)
      }
    }
  }

  return(output)
}




delta.method.acre <- function(cov.linked, param.values, link.funs = NULL, new.covariates = NULL,
                                 name.og = NULL, name.extend = NULL, par.ext.cov.provided = NULL, df.param = NULL,
                                 gam.model.full = NULL, gam.output = NULL, back.trans = TRUE){
  ## Delta method: G(x) is a vector of n functions, where x is a vector of m variables;
  ## then VAR(G(x)) = G'(x) %*% VAR(x) %*% t[G'(x)].
  ## VAR(x) is an m*m matrix, G'(x) is an n*m matrix, and G'(x)[i, j] = d(g[i])/d(x[j]).

  ## In this case, m is always the number of rows or columns of cov.linked, which equals the length of param.values.
  m <- length(param.values)

  ## In the simple case, no new covariates are involved, and "link.funs" should be provided;
  ## otherwise, "new.covariates", "name.og", "df.param" will be used together, and "link.funs" should be NULL.
  if(!is.null(link.funs)){
    ## In this case, G'(x) is an n*n diagonal matrix.
    n <- length(link.funs)
    G_grad <- matrix(0, nrow = n, ncol = n)
    for(i in 1:n){
      if(back.trans){
        link <- link.funs[i]
      } else {
        link <- 'identity'
      }

      x <- param.values[i]
      if(link == 'identity'){
        G_grad[i, i] <- 1
      } else if(link == 'log'){
        G_grad[i, i] <- exp(x)
      } else if(link == 'logit'){
        G_grad[i, i] <- exp(x)/(exp(x) + 1)^2
      }

    }

  } else {
    ## if(nrow(new.covariates) != 1) stop('Only 1 row in new.covariates is accepted.')
    ## If "new.covariates" is provided, then we are going to calculate the covariance matrix
    ## for all original parameters used in this model. Let "n" be length(unique(name.og));
    ## then G'(x) will be an n * m matrix.
    unique_name_og <- unique(name.og)
    n <- length(unique_name_og)
    G_grad <- vector('list', n)
    n_new_df <- nrow(new.covariates)
    ## G_grad = matrix(0, nrow = n, ncol = m)


    for(i in 1:n){
      ## For each row, the function is G_i(x) = h^(-1)(t(w)x), where w is the covariate values,
      ## x is the estimated parameter values, and h^(-1)(.) is the inverse of the link function.
      par_name <- unique_name_og[i]
      ## Extract the indices of this parameter in name.og, which holds the
      ## original parameter names of each column of cov.linked.
      index_p <- which(name.og == par_name)
      x <- param.values[index_p]
      len_p <- length(index_p)

      tem <- subset(df.param, df.param$par == par_name)
      if(back.trans){
        ## When a parameter is extended but not all covariates are provided, set the link to 'identity' as well,
        ## because in this case we will show, for example, "D_int_link" and "D_beta1_link".
        if((par_name %in% name.extend) & !(par_name %in% par.ext.cov.provided)){
          link <- 'identity'
        } else {
          link <- tem[, 'link']
        }

      } else {
        ## If we do not back-transform, it means, for example, that if "D" has "D_int" and "D_beta1",
        ## then we just use the delta method to calculate var(D_int + D_beta1 * x1), instead of
        ## var(exp(D_int + D_beta1 * x1)).
        link <- 'identity'
      }

      n_col_full <- tem[, 'n_col_full']
      n_col_mask <- tem[, 'n_col_mask']

      ## Let w be the covariate values for each parameter.
      ## w = numeric(n_col_full + n_col_mask)

      if((par_name %in% name.extend) & (par_name %in% par.ext.cov.provided)){
        G_grad[[i]] <- matrix(0, nrow = n_new_df, ncol = m)
        w <- matrix(0, nrow = n_new_df, ncol = n_col_full + n_col_mask)
        if(n_col_full > 1){
          gam.model <- gam.output[[par_name]][["gam_non_mask"]]
          w[, 1:n_col_full] <- get.DX.new.gam(gam.model, new.covariates)
        } else {
          w[, 1] <- matrix(1, nrow = n_new_df, ncol = 1)
        }

        if(n_col_mask > 0){
          gam.model <- gam.output[[par_name]][["gam_mask"]]
          tem <- get.DX.new.gam(gam.model, new.covariates)
          ## Get rid of the first column since the intercept is not here.
          w[, (n_col_full + 1):(n_col_full + n_col_mask)] <- tem[, -1]
        }

      } else {
        G_grad[[i]] <- matrix(0, nrow = len_p, ncol = m)
        w <- diag(len_p)
      }


      n_row <- nrow(w)

      if(link == 'identity'){
        for(k in 1:n_row){
          for(j in 1:len_p){
            G_grad[[i]][k, index_p[j]] <- w[k, j]
          }
        }
      } else if(link == 'log'){
        for(k in 1:n_row){
          common_component <- exp(sum(w[k,] * x))
          for(j in 1:len_p){
            G_grad[[i]][k, index_p[j]] <- common_component * w[k, j]
          }
        }

      } else if(link == 'logit'){
        for(k in 1:n_row){
          common_component <- exp(sum(w[k,] * x))
          for(j in 1:len_p){
            G_grad[[i]][k, index_p[j]] <- common_component * w[k, j] / (common_component + 1) ^ 2
          }
        }
      }

    }


    G_grad <- do.call('rbind', G_grad)
  }
  return(G_grad %*% cov.linked %*% t(G_grad))

}


## Since DX could be very large, we only calculate and return the diagonal elements.
delta.for.pred <- function(DX, par.value.linked, vcov.matrix, log.scale){
  ## Only used for D, so only consider the situation where the link function is log.
  n <- nrow(DX)
  p <- ncol(DX)

  if(!log.scale){
    G <- exp(DX %*% par.value.linked)
    G <- matrix(rep(G, p), ncol = p)
    G <- G * DX
  } else {
    G <- DX
  }


  output <- numeric(n)
  for(i in 1:n){
    tem_G <- G[i, , drop = FALSE]
    output[i] <- (tem_G %*% vcov.matrix %*% t(tem_G))[1,1]
  }

  return(output)

}


fixed.par.vcov.add <- function(m, p, type){
  ## When p is NULL, n_added will be zero.
  n_added <- length(p)
  if(n_added == 0){
    output <- m
  } else {
    if(type == 'linked') p <- paste(p, "link", sep = "_")

    m_new <- matrix(0, nrow = nrow(m) + n_added, ncol = ncol(m) + n_added)

    new_name <- c(rownames(m), p)
    rownames(m_new) <- new_name
    colnames(m_new) <- new_name
    m_new[1:nrow(m), 1:ncol(m)] <- m

    dim_m <- data.frame(par = ori.name(new_name), dim_name = new_name, order_og = 1:nrow(m_new))
    full_par <- fulllist.par.generator()
    full_par <- data.frame(par = full_par, o = 1:length(full_par))
    dim_m <- merge(dim_m, full_par, by = 'par')
    dim_m <- dim_m[order(dim_m$o, dim_m$order_og),,drop = FALSE]

    output <- matrix(0, nrow = nrow(m_new), ncol = ncol(m_new))
    rownames(output) <- dim_m$dim_name
    colnames(output) <- dim_m$dim_name

    tem <- matrix(0, nrow = nrow(m_new), ncol = ncol(m_new))
    for(i in 1:nrow(m_new)){
      i_m_new <- dim_m$order_og[i]
      tem[i, ] <- m_new[i_m_new, ]
    }

    for(i in 1:ncol(m_new)){
      i_m_new <- dim_m$order_og[i]
      output[, i] <- tem[, i_m_new]
    }

  }

  return(output)

}







vector.to.df <- function(vec){
  name <- names(vec)
  if(is.null(name)) name <- seq(length(vec))
  output <- data.frame(name = name, value = as.vector(vec))
  return(output)
}

gaussian.confint.cal <- function(object, types, pars, new.covariates, q.lower, q.upper){

  output <- vector('list', length(types))
  names(output) <- types
  for(i in types){
    ## if(i == 'linked'){
    ##  df_est = vector.to.df(coef.acre(object, types = 'linked', pars = pars))
    ##  df_std = vector.to.df(stdEr.acre(object, types = 'linked', pars = pars))
    ## }

    if(i == 'fitted' | i == 'linked'){
      ## The "fitted" interval is just back-transformed from the "linked" confidence interval.
      df_est <- vector.to.df(coef.acre(object, types = 'linked', pars = pars, new.covariates = new.covariates))
      df_std <- vector.to.df(stdEr.acre(object, types = 'linked', pars = pars, new.covariates = new.covariates, show.fixed.par = FALSE))
    }

    if(i == 'derived'){
      df_est <- vector.to.df(coef.acre(object, types = 'derived'))
      df_std <- vector.to.df(stdEr.acre(object, types = 'derived', show.fixed.par = FALSE))
    }
    colnames(df_est) <- c('par', 'est')
    df_est$std <- 0
    u_name <- unique(df_est$par)
    tem <- vector('list', length(u_name))
    names(tem) <- u_name
    for(k in u_name){
      tem_est <- subset(df_est, df_est$par == k)
      tem_std <- subset(df_std, df_std$name == k)
      if(nrow(tem_std) > 0) tem_est$std <- tem_std$value
      tem[[k]] <- tem_est
    }
    tem <- do.call('rbind', tem)
    tem$lower <- tem$est + q.lower * tem$std
    tem$upper <- tem$est + q.upper * tem$std


    if(i == 'fitted') tem$par <- gsub("_link", "", tem$par)

    output[[i]] <- tem

  }

  return(output)
}


ori.name <- function(char){
  ## This order matters: we must try to match 'sigma.b0.ss' and 'sigma.toa' first,
  ## otherwise "sigma.b0.ss.(Intercept)" will be recognised as 'sigma',
  ## so we must put 'sigma' after any 'sigma.xx' parameters.
  fulllist.par <- c('g0', 'lambda0', 'z', 'shape.1',
                   'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
                   'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa',
                   "sigma.b0.ss", 'D', 'sigma', 'mu')
  ## The pattern is the original parameter name at the start.
  fulllist_reg_exp <- gsub("\\.", "\\\\\\.", fulllist.par)
  pattern <- paste("^", fulllist_reg_exp, sep = "")
  n <- length(char)
  output <- character(n)
  for(i in 1:n){
    pattern_matched <- FALSE

    ## First try to match a pattern like "kappa.".
    for(j in 1:length(pattern)){
      tem <- regexpr(pattern[j], char[i])
      if(tem == 1){
        pattern_matched <- TRUE
        output[i] <- fulllist.par[j]
        break
      }
    }

    ## If we cannot find a pattern like "kappa.", then check whether the input is just the original parameter name.

    if(!pattern_matched){
      for(j in 1:length(fulllist.par)){
        if(char[i] == fulllist.par[j]){
          pattern_matched <- TRUE
          output[i] <- fulllist.par[j]
          break
        }
      }
    }

    if(!pattern_matched){
      msg <- paste0("Cannot find original parameter's name for the input: ", char[i])
      stop(msg)
    }

    ## End of loop for char[i].
  }

  return(output)
}

## Find the nearest mask point to an observed MRDS location.
find.nearest.mask <- function(locs, mask){
  dists <- distances(locs, mask)
  as.list(apply(dists, 1, function(x) which(x == min(x))[1]))
}

diy.mean <- function(vec){
  if(length(vec) > 0){
    if(is(vec, 'numeric')){
      return(mean(vec, na.rm = TRUE))
    } else {
      return(vec[1])
    }
  } else {
    if(is(vec, 'character')){
      return(character(0))
    } else {
      return(numeric(0))
    }
  }
}

## The weight methods are Shepard and Modified Shepard.
## The argument mask: the "mask" object generated by create.mask().

## The argument loc.cov: a list, data frame or matrix. If it is a list, each element should be a data frame or a matrix
## containing at least 3 columns, 2 of which are "x" and "y" for location, and the rest are the location-related
## covariates. A list allows different covariates from different sets of locations.

## The argument control.nn2: a list containing the arguments for the function RANN::nn2(); please refer to the help documentation in that
## package for further details.

## The argument control.weight: a list containing 3 elements that control the weighting method. The function RANN::nn2() gives us an M * k matrix,
## where M is the number of mask points, and each row contains k points provided by "loc.cov" which are nearest to this mask point.
## For these k nearest points, we use their values of each covariate to impute the estimated value at this mask point.
## Since the distances between this mask point and these k nearest points differ, we need a method to determine their weights.
## The argument "control.weight" is used to do so. Its 3 elements are "method", which could be "Shepard" or "Modified"; "q", which is used
## to control both methods; and "r", which is used to control the "Modified" method.
## The "Shepard" method: w ~ (1/d)^q.
## The "Modified" method: the modified Shepard method, w ~ [max(0, r - d)/(r * d)]^q.


## The argument control.char: a list containing 2 elements (weighted_nn, k), used to control character-type covariates.
## 1st element - "weighted_nn": a logical value, FALSE by default. If TRUE, use the value with the highest weight as the value for each mask point;
## for example, using RANN::nn2() and the control.weight argument, a mask point has 3 nearest points with w = c(0.4, 0.3, 0.3), and the values of these 3 points
## are ("a", "b", "b"); then "a" has a weight of 0.4 and "b" has 0.6 in total, so the value assigned to this mask point will be "b".
## If FALSE, then the KNN method will be used.

## 2nd element - k: used for the KNN method when weighted_nn is FALSE.

location.cov.to.mask <- function(mask, loc.cov, control.nn2 = NULL, control.weight = NULL, control.char = NULL){
  ## Declare the function nn2 from RANN.
  f <- RANN::nn2

  stopifnot(is(mask, 'list'))
  n.sessions <- length(mask)
  n.masks <- sapply(mask, nrow)

  if(is(loc.cov, 'data.frame') | is(loc.cov, 'matrix')){
    stopifnot(all(c('x', 'y') %in% colnames(loc.cov)))
    if("session" %in% colnames(loc.cov)){
      stopifnot(ncol(loc.cov) > 3)
      loc.cov <- loc.cov[!duplicated(loc.cov[,c('session', 'x', 'y')]), ]
      name_cov <- colnames(loc.cov)[-which(colnames(loc.cov) %in% c('session', 'x', 'y'))]
    } else {
      stopifnot(ncol(loc.cov) > 2)
      loc.cov <- loc.cov[!duplicated(loc.cov[,c('x', 'y')]), ]
      name_cov <- colnames(loc.cov)[-which(colnames(loc.cov) %in% c('x', 'y'))]
    }

    n_loc_cov <- 1
    ## This index indicates which covariate is extracted from which element of loc.cov.
    index_name_cov <- rep(1, length(name_cov))
    loc.cov <- list(loc.cov)
  } else {
    stopifnot(is(loc.cov, 'list'))
    n_loc_cov <- length(loc.cov)
    name_cov <- vector('list', n_loc_cov)
    ## This index indicates which covariate is extracted from which element of loc.cov.
    index_name_cov <- vector('list', n_loc_cov)
    for(i in 1:n_loc_cov){
      tem <- loc.cov[[i]]
      stopifnot(any(is(tem, 'data.frame'), is(tem, 'matrix')))
      stopifnot(all(c('x', 'y') %in% colnames(tem)))
      if('session' %in% colnames(tem)){
        stopifnot(ncol(tem) > 3)
        tem <- tem[!duplicated(tem[,c('session', 'x', 'y')]), ]
        name_cov[[i]] <- colnames(tem)[-which(colnames(tem) %in% c('session', 'x', 'y'))]
      } else {
        stopifnot(ncol(tem) > 2)
        tem <- tem[!duplicated(tem[,c('x', 'y')]), ]
        name_cov[[i]] <- colnames(tem)[-which(colnames(tem) %in% c('x', 'y'))]
      }
      index_name_cov[[i]] <- rep(i, length(name_cov[[i]]))
      loc.cov[[i]] <- tem
    }

    name_cov <- do.call('c', name_cov)
    index_name_cov <- do.call('c', index_name_cov)
  }


  if(any(duplicated(name_cov))) stop('duplicated covariates found.')


  ## Set default values
  ## when the user does not assign any value to these arguments.
  if(is.null(control.nn2)) control.nn2 <- list(k = 10)
  if(is.null(control.weight)) control.weight <- list(q = 2, method = 'Shepard', r = NULL)
  if(is.null(control.char)) control.char <- list(weighted_nn = FALSE, k = 1)

  ## When the user only assigns some components of these arguments.
  if(is.null(control.nn2$k)) control.nn2$k <- 10
  if(is.null(control.char$weighted_nn)) control.char$weighted_nn <- FALSE
  if(is.null(control.char$k)) control.char$k <- 1
  if(is.null(control.weight$method)) control.weight$method <- 'Shepard'
  if(is.null(control.weight$q)) control.weight$q <- 2

  output <- vector('list', n.sessions)

  for(s in 1:n.sessions){
    tem_mask <- mask[[s]]
    output[[s]] <- data.frame(session = s, mask = 1:n.masks[s])

    o_nn2 <- vector('list', n_loc_cov)
    w <- vector('list', n_loc_cov)

    loc_cov_session_subset <- vector('list', n_loc_cov)

    ## Each i denotes one provided loc.cov data frame, which might contain multiple covariates.
    for(i in 1:n_loc_cov){
      tem_cov <- loc.cov[[i]]
      if('session' %in% colnames(tem_cov)){
        tem_cov <- subset(tem_cov, tem_cov$session == s)
        if(nrow(tem_cov) == 0) stop('session has been provided in loc.cov, please make sure all sessions are included.')
        tem_cov <- tem_cov[,-which(colnames(tem_cov) == 'session')]
      }

      loc_cov_session_subset[[i]] <- tem_cov

      tem_nn2_par <- control.nn2
      tem_nn2_par$k <- min(tem_nn2_par$k, nrow(tem_cov))
      tem_nn2_par$data <- tem_cov[, c('x', 'y')]
      tem_nn2_par$query <- tem_mask

      ## Ideally, the range denoted by tem_nn2_par$data should cover all of tem_nn2_par$query, or
      ## at least the boundary of tem_nn2_par$query should not exceed tem_nn2_par$data by too much;
      ## here "too much" means 20% of the range of the xlim or ylim of the data.
      xlim_data <- range(tem_nn2_par$data[, 'x'])
      xscale_data <- diff(xlim_data)
      ## When the x coordinates of the data are all the same, e.g. the data are provided along a line,
      ## I'm not sure how to determine the "scale", so just take abs(x) for now,
      if(xscale_data == 0) xscale_data <- abs(xlim_data[1])
      xlim_data_max <- c(xlim_data[1] - 0.7 * xscale_data, xlim_data[2] + 0.7 * xscale_data)
      xlim_query <- range(tem_nn2_par$query[,'x'])

      ylim_data <- range(tem_nn2_par$data[, 'y'])
      yscale_data <- diff(ylim_data)
      if(yscale_data == 0) yscale_data <- abs(ylim_data[1])
      ylim_data_max <- c(ylim_data[1] - 0.7 * yscale_data, ylim_data[2] + 0.7 * yscale_data)
      ylim_query <- range(tem_nn2_par$query[,'y'])

      o_nn2[[i]] <- do.call('f', tem_nn2_par)

      ## Remove zero indices, which may occur when we do a "radius" search.

      ## The element nn.idx contains the indices in loc.cov of the k nearest points for each mask point;
      ## nn.dists contains the distances between each mask point and those k nearest points;
      ## both of them are M * k matrices, where M is the number of mask points.


      zero_index <- which(apply(o_nn2[[i]]$nn.idx, 2, sum) == 0)
      if(length(zero_index) > 0){
        o_nn2[[i]]$nn.idx <- o_nn2[[i]]$nn.idx[, -zero_index, drop = FALSE]
        o_nn2[[i]]$nn.dists <- o_nn2[[i]]$nn.dists[,-zero_index, drop = FALSE]
        if(ncol(o_nn2[[i]]$nn.idx) == 0){
          stop('Please specify a larger radius for radius searching.')
        }
      }


      ## For each mask point, the value of each covariate is a weighted average of all the provided covariate values;
      ## there are two ways to determine the weights:
      ## Shepard: w ~ (1/d)^q, where q is a provided constant;
      ## Modified Shepard: w ~ [max(0, r - d) / (r * d)] ^ q, where r and q are provided constants.

      ## Calculate the weight matrix (for any row of w, sum(w[i,]) === 1).
      if(control.weight$method == 'Shepard'){
        w[[i]] <- (1/o_nn2[[i]]$nn.dists) ^ control.weight$q
        w[[i]] <- w[[i]] / apply(w[[i]], 1, sum)
      } else if(control.weight$method == 'Modified'){
        if(is.null(control.weight$r)){
          stop('please specify a valid radius for Modified Shepard Weight by control.weight$r = xx')
        }
        w[[i]] <- (max(0, control.weight$r - o_nn2[[i]]$nn.dists) / (control.weight$r * o_nn2[[i]]$nn.dists)) ^ control.weight$q
        w[[i]] <- w[[i]] / apply(w[[i]], 1, sum)
      } else {
        stop('weight method only suppors "Shepard" and "Modified" currently.')
      }

      ## If any location with covariates lies exactly on a mask point, the weight is NaN, so convert it to 1.
      if(any(is.nan(w[[i]]))){
        w[[i]][is.nan(w[[i]])] <- 1
      }
      ## End of for(i in 1:n_lco_cov).
    }




    for(j in name_cov){
      ## First, locate the index of the element we got this covariate from.
      i <- index_name_cov[which(j == name_cov)]
      values <- loc_cov_session_subset[[i]][, j]

      if(is.numeric(values)){
        ## Generate a matrix with the same dimension as o_nn2[[i]]$nn.idx, and
        ## mat_values[i, j] = values[o_nn2[[i]]$nn.idx[i, j]]
        mat_values <- matrix(values[o_nn2[[i]]$nn.idx], nrow = n.masks[s])
        output[[s]][[j]] <- apply(mat_values * w[[i]], 1, sum)
      } else {
        values <- as.character(values)
        ## Generate a matrix with the same dimension as o_nn2[[i]]$nn.idx, and
        ## mat_values[i, j] = values[o_nn2[[i]]$nn.idx[i, j]]
        mat_values <- matrix(values[o_nn2[[i]]$nn.idx], nrow = n.masks[s])
        v <- character(n.masks[s])

        if(control.char$weighted_nn){
          ## For each row, select the element from the 'nearest neighbours' with the largest weight.
          for(n in 1:n.masks[s]){
            tem_v <- mat_values[n,]
            tem_w <- w[[i]][n,]
            tem_agg <- aggregate(tem_w, list(v = tem_v), sum)
            v[n] <- tem_agg$v[which.max(tem_agg$x)]
          }
        } else {
          ## If not using the weighted closest point, we can simply use the KNN method for classification.

          ## Since KNN is based on nn2's result, the k here must be smaller than or equal to the number of
          ## "closest points" generated by nn2.
          stopifnot(control.char$k <= ncol(o_nn2[[i]]$nn.idx))

          ## Drop the (k + 1)th and subsequent columns.
          mat_values <- mat_values[, 1:control.char$k, drop = FALSE]


          if(control.char$k != 1){
            ## For each row, select the most frequent value among its k closest neighbours.
            v <- apply(mat_values, 1, function(x) names(which.max(table(x))))
          } else {
            ## If k == 1, then we only take the nearest neighbour, and mat_values is an n*1 matrix, which is easy.
            v <- as.vector(mat_values)


          }
        }

        output[[s]][[j]] <- v
      }
      ## End of covariate j.
    }

    ## End of session s.
  }
  return(do.call('rbind', output))

}

##
distance.cov.to.mask <- function(mask, dist.cov) {
  ## Directly calculate distance covariates for the provided mask.
  if(!is.list(dist.cov)){
    stop("dist.cov must be a list containing distance covariates")
  }

  ## Directly calculate the distance from each mask point to the nearest covariate
  ## area.
  calculated_dist.cov <- convert.dist.cov.to.loc.cov(dist.cov, NULL, mask)

  return(calculated_dist.cov)
}

par.extend.create <- function(loc.cov = NULL, mask = NULL, convert.loc2mask = list(),
                             dist.cov = NULL, session.cov = NULL, trap.cov = NULL, time.loc.cov = NULL){


  if(any(!is.null(loc.cov), !is.null(session.cov), !is.null(trap.cov), !is.null(dist.cov), !is.null(time.loc.cov))){
    par.extend <- list()

    if(!is.null(session.cov)) stopifnot(is(session.cov, 'data.frame'))
    if(!is.null(trap.cov)) stopifnot(is(trap.cov, 'data.frame'))

    ## For covariates with an area edge, the nearest distance to each mask point will be assigned
    ## as a mask-level covariate. Here we calculate the distances and assign them to loc.cov; then
    ## they will be handled together with the other location-related covariates.

    loc.cov <- convert.dist.cov.to.loc.cov(dist.cov = dist.cov, loc.cov = loc.cov, mask = mask)

    loc.cov <- convert.time.loc.cov.to.loc.cov(time.loc.cov = time.loc.cov, loc.cov = loc.cov, session.cov = session.cov)
    ## If location-related covariates are provided, convert them to a mask-level data frame.
    if(!is.null(loc.cov)){
      convert.loc2mask$loc.cov <- loc.cov
      convert.loc2mask$mask <- mask

      mask_cov <- do.call('location.cov.to.mask', convert.loc2mask)
    } else {
      mask_cov <- NULL
    }



    par.extend$data <- list()
    par.extend$data$session <- session.cov
    par.extend$data$trap <- trap.cov
    par.extend$data$mask <- mask_cov
  } else (
    par.extend <- NULL
  )

  return(list(output = par.extend, loc.cov = loc.cov))
}


demo.traps <- function(){
  o <- data.frame(x = rep(c(0,5) ,each = 3), y = rep(c(0, 5, 10), 2))
  return(o)
}

demo.traps2 <- function(){
  o <- data.frame(x = rep(20, 3), y = c(20, 25, 30))
  return(o)
}

demo.session.cov <- function(){
  o <- data.frame(session = 1:3, weather = c('sunny', 'rainy', 'sunny'))
  return(o)
}

demo.trap.cov <- function(){
  o <- data.frame(trap = 1:6, brand = rep(c('sony', 'panasonic'), each = 3))
  return(o)
}

demo.trap.cov2 <- function(){
  o <- data.frame(trap = 1:3, brand = c('sony', 'panasonic', 'panasonic'))
  return(o)
}

demo.loc.cov <- function(){
  o <- data.frame(x = rep(c(-20, 2.5, 25), each = 3), y = rep(c(-20, 5, 30), 3),
                 noise = c(6, 10, 11, 7, 12, 10, 11, 9, 8),
                 forest_volumn = c(rep('high', 2), rep('median', 3), rep('low', 4)))
  return(o)
}


sim.args.generator <- function(sim.name){
  if (!sim.name %in% get.dataset.names()) {
    stop('invalid input for "sim.name", which should be one of the following: "',
         paste(get.dataset.names(), collapse = '", "'), '"')
  }

  output <- list()
  ## Generate some common settings.
  traps <- demo.traps()
  control.mask <- list(buffer = 30)
  model <- NULL
  survey.length <- NULL
  ss.opts <- NULL
  cue.rates <- NULL
  n.sessions <- NULL
  session.cov <- demo.session.cov()
  trap.cov <- demo.trap.cov()
  loc.cov <- demo.loc.cov()


  if(sim.name == 'dist_hn'){
    param <- list(g0 = 1, sigma = 5.33, alpha = 5.62, D = 2300)
    detfn <- 'hn'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
  } else if(sim.name == 'bearing_dist_hn'){
    param <- list(g0 = 1, sigma = 5.39, kappa = 51.12, alpha = 5.02, D = 2265)
    detfn <- 'hn'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
  } else if(sim.name == 'bearing_hn'){
    param <- list(g0 = 1, sigma = 5.53, kappa = 53.88, D = 2172)
    detfn <- 'hn'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
  } else if(sim.name == 'ihd'){
    param <- list(g0 = 1, sigma = 5.38, D = c(8.3, -0.2, -0.1))
    detfn <- 'hn'
    model <- list(D = ~forest_volumn)
    session.cov <- NULL
    trap.cov <- NULL
  } else if(sim.name == 'ihd_ext'){
    param <- list(g0 = 1, sigma = c(1.95, -0.25), D = c(8.3, -0.1))
    detfn <- 'hn'
    model <- list(D = ~noise, sigma = ~brand)
    session.cov <- NULL
  } else if(sim.name == 'mul_ses'){
    param <- list(g0 = 0.7135, sigma = 3.3, kappa = 14.8, alpha = 3.77, D = 2533)
    traps <- list(traps, demo.traps2())
    control.mask <- list(buffer = 15)
    detfn <- 'hn'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL

  } else if(sim.name == 'mul_ses_ext'){
    param <- list(g0 = c(1.172, -0.8), sigma = c(1.1383, 0.1407), kappa = 14.8, alpha = 3.77, D = 2533)
    traps <- list(traps, demo.traps2())
    control.mask <- list(buffer = 15)
    detfn <- 'hn'
    session.cov <- session.cov[1:2,]
    trap.cov <- rbind(trap.cov, demo.trap.cov2())
    trap.cov$session <- c(rep(1, 6), rep(2, 3))
    loc.cov <- NULL
    model <- list(g0 = ~weather, sigma = ~brand)

  } else if(sim.name == 'simple_hhn'){
    param <- list(sigma = 3.66, lambda0 = 4.29, D = 2657)
    detfn <- 'hhn'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL

  } else if(sim.name == 'hhn_cue'){
    param <- list(sigma = 3.66, lambda0 = 4.29, D = 500)
    detfn <- 'hhn'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
    cue.rates <- c(5, 4, 11, 2, 3)
    survey.length <- 1
  } else if(sim.name == 'simple_hr'){
    param <- list(g0 = 0.88, sigma = 7, z = 7.5, D = 2665)
    detfn <- 'hr'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
  } else if(sim.name == 'ss'){
    param <- list(b0.ss = 89.07, b1.ss = 4.13, sigma.ss = 9.52, D = 2657)
    detfn <- 'ss'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
    ss.opts <- list(cutoff = 60)
  } else if(sim.name == 'ss_toa'){
    param <- list(b0.ss = 90, b1.ss = 4, sigma.ss = 9.34, sigma.toa = 0.002, D = 2430)
    detfn <- 'ss'
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
    ss.opts <- list(cutoff = 60)
  } else if(sim.name == 'ind_bearing_dist'){
    param <- list(g0 = 0.68, sigma = 9.8, kappa = 10, alpha = c(0.56, 0.16), D = c(4.15, 0.62), mu = 8.8)
    detfn <- 'hn'
    control.mask <- list(buffer = 15)
    loc.cov <- NULL
    model <- list(D = ~weather, alpha = ~brand)
    n.sessions <- 3
    survey.length <- c(1, 2, 1)
  } else if(sim.name == 'ind_toa_hhn'){
    param <- list(sigma = 2.13, lambda0 = 8.33, sigma.toa = 0.0011, D = c(6.81, -0.12), mu = 8.2)
    detfn <- 'hhn'
    control.mask <- list(buffer = 15)
    session.cov <- NULL
    trap.cov <- NULL
    model <- list(D = ~noise)
    n.sessions <- 2

  } else if(sim.name == 'ind_ss'){
    param <- list(b0.ss = c(4.3, 0.1), b1.ss = 0.8872, sigma.ss = 10.1, D = c(4.3, -0.1), mu = 8.3)
    detfn <- 'ss'
    control.mask <- list(buffer = 15)
    session.cov <- NULL
    model <- list(D=~noise, b0.ss=~brand)
    n.sessions <- 3
    ss.opts <- list(cutoff = 60)
  } else if(sim.name == 'ind_ss_log'){
    ##################################################################################################
    ## The estimate seems to be biased; check it further.
    param <- list(b0.ss = 5.35, b1.ss = 0.08, sigma.ss = 9.788, D = 67, mu = 8.3)
    detfn <- 'ss'
    control.mask <- list(buffer = 15)
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
    n.sessions <- 3
    ss.opts <- list(cutoff = 75, ss.link = 'log')
    ####################################################################################################
  } else if(sim.name == 'ind_ss_sp'){
    ##################################################################################################
    ## The estimate seems to be biased; check it further.
    param <- list(b0.ss = 100.25, b1.ss = 0.082, sigma.ss = 9.788, D = 52, mu = 9.1)
    detfn <- 'ss'
    control.mask <- list(buffer = 15)
    session.cov <- NULL
    trap.cov <- NULL
    loc.cov <- NULL
    n.sessions <- 3
    ss.opts <- list(cutoff = 75, ss.link = 'spherical')
    ##################################################################################################
  }


  ## Include all arguments that may differ from the defaults.
  output$detfn <- detfn
  output$param <- param
  output$model <- model
  output$traps <- traps
  output$control.mask <- control.mask
  output$session.cov <- session.cov
  output$trap.cov <- trap.cov
  output$loc.cov <- loc.cov
  output$survey.length <- survey.length
  output$ss.opts <- ss.opts
  output$cue.rates <- cue.rates
  output$n.sessions <- n.sessions

  return(output)

}


fit.args.generator.from.sim <- function(sim.name, fit.args){
  output <- fit.args

  if(sim.name == 'dist_hn'){
    output$fix <- list(g0 = 1)
  } else if(sim.name == 'bearing_dist_hn'){
    output$fix <- list(g0 = 1)
  } else if(sim.name == 'bearing_hn'){
    output$fix <- list(g0 = 1)
  } else if(sim.name == 'ihd'){
    output$fix <- list(g0 = 1)
  } else if(sim.name == 'ihd_ext'){
    output$fix <- list(g0 = 1)
  } else if(sim.name == 'mul_ses'){
    output$sv <- list(kappa = 20)
  } else if(sim.name == 'mul_ses_ext'){
    output$sv <- list(kappa = 20)
  } else if(sim.name == 'simple_hhn'){
    ## Nothing needs to be done.
  } else if(sim.name == 'hhn_cue'){
    ## Nothing needs to be done.
  } else if(sim.name == 'simple_hr'){
    output$sv <- list(z = 5)
  } else if(sim.name == 'ss'){
    output$sv <- list(b0.ss = 90, b1.ss = 4, sigma.ss = 10)
  } else if(sim.name == 'ss_toa'){
    output$sv <- list(b0.ss = 90, b1.ss = 4, sigma.ss = 10)
  } else if(sim.name == 'ind_bearing_dist'){
    ## Nothing needs to be done.
  } else if(sim.name == 'ind_toa_hhn'){
    ## Nothing needs to be done.
  } else if(sim.name == 'ind_ss'){
    ## Nothing needs to be done.
  } else if(sim.name == 'ind_ss_log'){
    ## Nothing needs to be done.
  } else if(sim.name == 'ind_ss_sp'){
    ## Nothing needs to be done.
  }

  return(output)
}


param.transform <- function(param, df.link, back = FALSE){
  par_names <- names(param)
  which_trans <- par_names[which(sapply(param, function(x) length(x) == 1))]

  if(back){
    for(i in which_trans){
      link <- df.link[which(df.link$par == i), 'link']
      param[[i]] <- unlink.fun(link, param[[i]])
    }
  } else {
    for(i in which_trans){
      link <- df.link[which(df.link$par == i), 'link']
      param[[i]] <- link.fun(link, param[[i]])
    }
  }


  return(param)
}


default.df.link <- function(){
  pars <- fulllist.par.generator()
  output <- data.frame(par = pars, link = character(length(pars)))
  for(i in 1:length(pars)){
    output[i,'link'] <- default.link(pars[i])
  }
  return(output)
}


diag.block.combine <- function(lst){
  n_lst <- sapply(lst, nrow)
  n <- sum(n_lst)
  out <- matrix(0, nrow = n, ncol = n)

  start_point <- 1

  for(i in 1:length(lst)){
    end_point <- start_point + n_lst[i] - 1
    indices <- seq(from = start_point, to = end_point)
    out[indices, indices] <- lst[[i]]
    start_point <- start_point + n_lst[i]
  }

  return(out)

}


convert.scale <- function(from, to){
  max_to <- max(to)
  min_to <- min(to)

  if(length(from) == 1){
    output <- 0.5 * (max_to + min_to)
  } else {
    min_from <- min(from)
    dist_og <- max(from) - min_from
    dist_to <- max_to - min_to
    if(dist_og == 0){
      output <- 0.5 * (max_to + min_to)
    } else {
      output <- (from - min_from) * dist_to / dist_og + min_to
    }
  }


  return(output)
}

circle.fun <- function(centre = data.frame(x = 0, y = 0), r = 0.5, npoints = 50){
  n <- nrow(centre)
  if(length(r) != n) r <- rep(r, length.out = n)
  output <- vector('list', n)
  for(i in 1:n){
    tt <- seq(0, 2*pi, length.out = npoints)
    xx <- centre[i, 1] + r[i] * cos(tt)
    yy <- centre[i, 2] + r[i] * sin(tt)
    output[[i]] <- data.frame(x = xx, y = yy, cir_index = i)
  }

  output <- do.call('rbind', output)

  return(output)
}


homo.digit <- function(x){
  x <- as.character(x)
  n_digit <- nchar(x)
  max_digit <- max(n_digit)

  n_zero <- max_digit - n_digit

  y <- character(length(x))
  for(i in 1:length(x)){
    y[i] <- paste(rep('0', n_zero[i]), collapse = "")
  }
  output <- paste(y, x, sep = "")
  return(output)
}




## Predict density using location coordinates; currently only used in the plot.

D.predict.for.plot <- function(fit, session.select = 1, new.data = NULL,
                              D.cov = NULL, xlim = NULL, ylim = NULL,
                              x.pixels = 50, y.pixels = 50, se.fit = FALSE, log.scale = FALSE, set.zero = NULL,
                              convert.loc2mask = NULL,
                              control.boot = list(correct.bias = FALSE, from.boot = TRUE), ...){
  if(!is.null(set.zero) & se.fit){
    warning('During the calculation of standard error, std_zero will be ignored.')
  }

  original_mask <- FALSE

  if(is.null(new.data)){

    ## If new.data is not provided, but xlim and ylim are, we use xlim and ylim to build the mask instead of
    ## extracting the mask from the fit.
    if(any(!is.null(xlim), !is.null(ylim))){
      if(!all(!is.null(xlim), !is.null(ylim))) stop('please provide both xlim and ylim.')
      x <- seq(from = xlim[1], to = xlim[2], length.out = x.pixels)
      y <- seq(from = ylim[1], to = ylim[2], length.out = y.pixels)
      mask <- data.frame(x = rep(x, each = y.pixels), y = rep(y, x.pixels))
    } else {
      mask <- get.mask(fit)[[session.select]]
      original_mask <- TRUE
    }

  } else {
    ## If new.data is provided, then xlim and ylim just do what they are supposed to do,
    ## i.e. trim the plot, instead of building a "mask".
    stopifnot(any(is(new.data, 'data.frame'), is(new.data, 'matrix')))
    stopifnot(all(c('x', 'y') %in% colnames(new.data)))
    mask <- new.data[, c('x', 'y'), drop = FALSE]
  }


  if('D' %in% get.par.extend.name(fit)){

    ## Build old_covariates;
    ## we interpolate it again whether or not there is a new mask grid because, in theory, the user
    ## could use the same mask grid but a different convert.loc2mask.

    old_covariates <- as.data.frame(mask)

    if(original_mask&is.null(convert.loc2mask)&is.null(D.cov)){
      ## When there is no sign of any modification, set old_loc_cov to NULL to make sure
      ## the function enters the process that tries to obtain mask-level data, to avoid re-interpolation.
      old_loc_cov <- NULL
    } else {
      old_loc_cov <- get.loc.cov(fit)
      ## When we cannot find loc.cov in the model fitting object, it is possible that D is
      ## not extended at the mask level; it is also possible that the model fitting object comes
      ## from a simulation study.
      mask_level_dat_extract <- FALSE
    }

    ## If any of the following holds:
    ##
    ## - the fit does not contain any location covariate information,
    ## - the mask we will be predicting over has changed from the original fit,
    ## - the covariate information we will be predicting with has changed from
    ## the original fit.
    ##
    ## then we will re-interpolate the location covariates.
    if(is.null(old_loc_cov)){

      ## Grab the original parameter interpolation data.
      old_loc_cov <- get.par.extend.data(fit)$mask

      ## If NULL, then the model was fitted without location covariates.
      if(!is.null(old_loc_cov)){

        mask_level_dat_extract <- TRUE

        ## Grab the correct session.
        if('session' %in% colnames(old_loc_cov)) {
          old_loc_cov <- subset(old_loc_cov,
                               old_loc_cov$session == session.select)
        }

        ## When the locations/mask for prediction are given by `xlim`/`ylim` or
        ## `new.data`, we use the original mask-level data as `loc.cov`;
        ## otherwise, we do not need to do any conversion in the first place.
        if(!original_mask){
          if('session' %in% colnames(old_loc_cov)) {
            old_loc_cov <- old_loc_cov[, which(colnames(old_loc_cov)!='session'),
                                      drop = FALSE]
          }
          if('mask' %in% colnames(old_loc_cov)) {
            old_loc_cov <- old_loc_cov[, which(colnames(old_loc_cov)!='mask'),
                                      drop = FALSE]
          }
          tem_mask <- get.mask(fit)
          if(is(tem_mask, 'list')) tem_mask <- tem_mask[[session.select]]
          stopifnot(nrow(old_loc_cov) == nrow(tem_mask))
          old_loc_cov <- cbind(tem_mask, old_loc_cov)
        }
      }
    }

    if(!is.null(old_loc_cov)) {
      ## As mentioned just above, if we are using the original mask data and
      ## mask-level covariate data, we do not need to do any conversion at all.
      if(mask_level_dat_extract & original_mask) {
        cov_mask <- old_loc_cov
      } else {
        if(is.null(convert.loc2mask)){
          convert.loc2mask <- vector('list', 2)
          names(convert.loc2mask) <- c('mask', 'loc.cov')
        }
        ## We need to create this object regardless of whether convert.loc2mask
        ## is provided.
        convert.dist2mask <- vector('list', 2)
        names(convert.dist2mask) <- c('mask', 'dist.cov')

        convert.dist2mask$mask <- list(mask)
        convert.loc2mask$mask <- list(mask)

        ## We only need to interpolate location covariates, as we can directly
        ## calculate the value at each mask point for distance covariates.

        ## Separate distance covariates from location covariates.
        dist_cov_names <- names(fit$arg_input$dist.cov)
        separated_cov <- separate.dist.loc.cov(old_loc_cov, dist_cov_names)

        ## NOTICE: we aren't using separated_cov$distance; we want to calculate
        ## the value for each mask point directly, so we just want the distance covariate values
        ## contained within the fit object.
        convert.dist2mask$dist.cov <- fit$arg_input$dist.cov
        convert.loc2mask$loc.cov <- separated_cov$location

        loc_cov_mask <- do.call('location.cov.to.mask', convert.loc2mask)
        dist_cov_mask <- do.call('distance.cov.to.mask', convert.dist2mask)

        ## Combine the distance and location covariates into one data frame.
        ## The distance mask should (read: will) have x, y columns matching old_covariates.
        cov_mask <- cbind(loc_cov_mask, dist_cov_mask)
      }
      ## Make sure to remove any unnecessary session or mask columns.
      if(any(colnames(cov_mask) %in% c('session', 'mask'))){
        old_covariates <- cbind(old_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask', 'x', 'y')), drop = FALSE])
      } else {
        old_covariates <- cbind(old_covariates, cov_mask)
      }
    }

    session_data_in_model <- get.par.extend.data(fit)$session
    ## Session-related covariates are simple: just take them
    ## from the input arguments of the fit. But don't do this if new data
    ## are provided.
    if(!is.null(session_data_in_model) & is.null(new.data)){
      tem <- session_data_in_model[which(session_data_in_model$session == session.select), , drop = FALSE]
      tem <- tem[, -which(colnames(tem) == 'session'), drop = FALSE]
      for(i in colnames(tem)) old_covariates[[i]] <- tem[1,i]
    }

    if(!is.null(D.cov) | !is.null(new.data)){
      ## First deal with all scenarios in which new covariates may be provided.
      if(!is.null(D.cov)){
        stopifnot(is(D.cov, 'list'))
        stopifnot(any(c('session', 'location') %in% names(D.cov)))
      }

      ## Consider the possibility that the user may only want to change some covariates and keep the others the same
      ## as in the model; for example, with 3 covariates related to D, 'weather' (session-related), and 'noise' and 'forest_type' (location-related),
      ## the user may only want to change weather, or noise, or noise and forest_type, etc.

      ## So we create 2 data frames, one for all new covariates and one for old covariates. We must make sure these two data frames
      ## contain the same mask points, and then replace the covariates in the "old" data frame with the same columns in the "new" one,
      ## if the same covariate appears in the "new" data frame.
      if(!is.null(new.data)){
        ## In new.data, the user can include any location-related covariates directly, and it also contains x and y,
        ## so we can use it directly instead of the mask.
        new.covariates <- as.data.frame(new.data)
      } else {
        new.covariates <- as.data.frame(mask)
      }
      ## Build new.covariates based on all the information we have.
      if(!is.null(D.cov$location)){
        if(is.null(convert.loc2mask)){
          convert.loc2mask <- vector('list', 2)
          names(convert.loc2mask) <- c('mask', 'loc.cov')
        }
        convert.loc2mask$mask <- list(mask)
        convert.loc2mask$loc.cov <- D.cov$location

        cov_mask <- do.call('location.cov.to.mask', convert.loc2mask)
        new.covariates <- cbind(new.covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
      }


      ## Since we only plot one session, the number of rows of D.cov$session should be 1.
      if(!is.null(D.cov$session)){
        stopifnot(nrow(D.cov$session) == 1)
        for(i in colnames(D.cov$session)) new.covariates[[i]] <- D.cov$session[1,i]

      }
      ## Update old_covariates with new.covariates.

      for(i in colnames(new.covariates)){
        if(all(i != 'x', i != 'y')){
          old_covariates[[i]] <- new.covariates[[i]]
        }
      }
    }

    par_info <- get.data.param(fit)
    par_info <- subset(par_info, par_info$par == 'D')
    gam.model <- get.gam(fit, 'D')
    if(is(fit, 'acreboot')){
      values_link <- as.vector(coef(fit, types = 'linked', pars = 'D', correct.bias = control.boot$correct.bias))
    } else {
      values_link <- as.vector(coef(fit, types = 'linked', pars = 'D'))
    }

    values_link_og <- values_link

    if(!is.null(set.zero)) values_link[set.zero] <- 0

    tem <- get.extended.par.value(gam.model, par_info$n_col_full, par_info$n_col_mask,
                                 values_link, old_covariates, DX.output = TRUE)
    D.mask <- unlink.fun(link = par_info$link, value = tem$output)

    if(se.fit){
      DX <- tem$DX

      if(is(fit, 'acreboot')) {
        if(log.scale){
          type <- 'linked'
        } else {
          type <- 'fitted'
        }

        D.se <- as.vector(stdEr(fit, types = type, new.covariates = old_covariates, pars = 'D',
                               from.boot = control.boot$from.boot, show.fixed.par = FALSE))

      } else {
        vcov_matrix <- vcov(fit, types = 'linked', pars = 'D', show.fixed.par = FALSE)
        D.se <- sqrt(delta.for.pred(DX, values_link_og, vcov_matrix, log.scale))
      }


    }

  } else {
    ## When D is not extended, it is literally a constant, so just take the 'fitted' estimate of D from coef() and stdEr().
    D.mask <- rep(coef(fit, types = "fitted", pars = 'D'), nrow(mask))
    if(se.fit){
      if(log.scale){
        type <- "linked"
      } else {
        type <- "fitted"
      }

      D.se <- as.vector(rep(stdEr(fit, types = type, pars = 'D', show.fixed.par = FALSE), nrow(mask)))
    }

  }

  output <- as.data.frame(mask)
  output$est <- D.mask
  if(se.fit) output$std <- D.se


  return(output)
}



res.split <- function(res, n.sessions){
  n_pars <- ncol(res) - 1 - n.sessions
  ## Extract and remove the maximum gradient component.
  maxgrads <- res[, n_pars + 1]
  res_esa <- res[, (n_pars + 2):ncol(res), drop = FALSE]
  res <- res[, 1:n_pars, drop = FALSE]


  return(list(res = res, maxgrads = maxgrads, res_esa = res_esa))
}


var.from.res <- function(res, fixed = NULL){
  par_names <- colnames(res)

  if(!grepl('^esa', par_names[1])){
    to_be_rm <- which(ori.name(par_names) %in% fixed)
    if(length(to_be_rm) != 0){
      res <- res[,-to_be_rm, drop = FALSE]
      par_names <- par_names[-to_be_rm]
    }
  }

  n_pars <- ncol(res)

  se <- apply(res, 2, "sd", na.rm = TRUE)
  names(se) <- par_names

  corr <- diag(n_pars)
  dimnames(corr) <- list(par_names, par_names)
  vcov <- diag(n_pars)
  diag(vcov) <- se^2
  dimnames(vcov) <- list(par_names, par_names)
  if(n_pars > 1){
    for (i in 1:(n_pars - 1)){
      for (j in (i + 1):n_pars){
        corr[i, j] <- corr[j, i] <- cor(res[, i], res[, j], use = "na.or.complete")
        vcov[i, j] <- vcov[j, i] <- corr[i, j]*se[i]*se[j]
      }
    }
  }

  return(vcov)

}


mce.from.res <- function(res, coefs, M, seed.mce){
  n_pars <- ncol(res)
  par_names <- colnames(res)
  bias <- apply(res, 2, "mean", na.rm = TRUE) - as.vector(coefs)
  ## Bootstrap to calculate the MCE for the bias and standard errors.
  bias.mce <- numeric(n_pars)
  se.mce <- numeric(n_pars)
  names(bias.mce) <- par_names
  names(se.mce) <- par_names
  if (M > 0){
    set.seed(seed.mce)
    converged <- which(!is.na(res[, 1]))
    n.converged <- length(converged)
    mce.boot <- matrix(sample(converged, size = n.converged*M,
                              replace = TRUE), nrow = M,
                       ncol = n.converged)
    for (i in par_names){
      par.boot <- matrix(res[mce.boot, i], nrow = M, ncol = n.converged)
      bias.mce[i] <- sd(apply(par.boot, 1, mean) - coefs[i] - bias[i])
      se.mce[i] <- sd(apply(par.boot, 1, sd))
    }
  } else {
    bias.mce <- NA
    se.mce <- NA
  }

  return(list(bias = bias, bias.mce = bias.mce, se.mce = se.mce))
}


## Used to resolve the determination of "types" and "pars" in various methods.
types.pars.sol <- function(types, pars, new.covariates){

  ## If new.covariates is provided, we need "fitted".
  ## if(!is.null(new.covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  ## If types is still empty, set it to 'linked'.
  if(is.null(types)) types <- 'linked'

  ## Confirm that types is within the valid set of options.
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))){
    stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  }

  if('esa' %in% pars){
    pars <- pars[-which(pars == 'esa')]
    ## If pars = 'esa' only, then treat it as NULL, and force types to be 'derived' only, whatever
    ## was assigned.
    if(length(pars) == 0){
      pars <- NULL
      types <- 'derived'
    } else {

      ## If the user assigns any other "pars" and only assigns types = "derived", add "linked" to types.
      if(length(types) == 1){
        if(types == 'derived') types <- c(types, 'linked')
      }
      if(!'derived' %in% types) types <- c(types, 'derived')
    }
  }

  ## After the 'esa' modification, if we still have assigned 'pars', but types is assigned as
  ## "derived" only, we need to add the default setting "linked" to it.
  ## Together with the 'esa' modification above, the basic logic is that "pars" has priority
  ## over "types".
  if(length(pars) > 0 & all(types == 'derived')) types <- c(types, 'linked')

  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }

  return(list(types = types, pars = pars))

}

################################################################################################################################


boot.res.to.CI <- function(res, level){
  n_par <- ncol(res)

  p_upper <- 0.5 + 0.5 * level
  p_lower <- 0.5 - 0.5 * level

  output <- matrix(0, nrow = n_par, ncol = 2)
  colnames(output) <- paste(c(p_lower, p_upper) * 100, "%", sep = " ")
  rownames(output) <- colnames(res)

  for(i in 1:n_par){
    output[i,] <- quantile(res[,i], c(p_lower, p_upper))
  }

  return(output)
}

## Transform the linked bootstrap results to the fitted type.
res.transform <- function(res, new.covariates, pars, object, back.trans = TRUE){
  name_extend <- get.par.extend.name(object)
  df_param <- get.data.param(object)
  col_name_og <- ori.name(colnames(res))
  output <- vector('list', length(pars))
  names(output) <- pars

  extend_par_covariates <- get.par.extend.covariate(object)

  if(is.null(name_extend) & !is.null(new.covariates)){
    warning('No parameter is extended, argument "new.covariates" will be ignored.')
    new.covariates <- NULL
  }

  name_extned_covariate_provided <- ext.par.in.new.df(name_extend, new.covariates, extend_par_covariates)


  for(j in pars){
    values_link <- res[, which(col_name_og == j), drop = FALSE]
    par_info <- subset(df_param, df_param$par == j)
    link <- par_info$link

    if(j %in% name_extend){

      if(!is.null(new.covariates) & (j %in% name_extned_covariate_provided)){
        gam <- get.gam(object, j)
        values_fitted <- get.extended.par.value(gam, par_info$n_col_full,
                                               par_info$n_col_mask, values_link, new.covariates,
                                               matrix.par.value = TRUE)
        colnames(values_fitted) <- rep(j, ncol(values_fitted))
        if(!back.trans){
          colnames(values_fitted) <- paste(colnames(values_fitted), "link", sep = "_")
        }
      } else {
          ## If no new covariates are assigned to this extended parameter (maybe there are no new.covariates at all,
          ## or at least one of the covariates for this parameter is not provided),
          ## then treat all of the relevant betas as identity-linked parameters.
          values_fitted <- values_link
          if(back.trans){
            colnames(values_fitted) <- gsub("_link$", "", colnames(values_fitted))
          }
          link <- 'identity'
      }


    } else {
      values_fitted <- values_link
      if(back.trans){
        colnames(values_fitted) <- j
      } else {
        colnames(values_fitted) <- paste(j, "link", sep = "_")
      }

    }

    if(back.trans){
      output[[j]] <- unlink.fun(link = link, value = values_fitted)
    } else {
      output[[j]] <- values_fitted
    }

  }

  output <- do.call('cbind', output)
  return(output)
}

## A method to do bias correction for the parametric bootstrap;
## see https://math.mit.edu/~dav/05.dir/class25-slides-all.pdf.
res.mod.for.CI <- function(res, est, correct.bias){
  if(correct.bias){
    res <- -1 * res
    res <- sweep(res, 2, 2 * est, FUN = "+")
  }
  return(res)
}



erf <- function(x){
  2*pnorm(x*sqrt(2)) - 1
}


convert.dist.cov.to.loc.cov <- function(dist.cov, loc.cov, mask){


  if(!is.null(dist.cov)){
    if(!is(dist.cov, 'list')) stop('dist.cov must be a named list.')
    n <- length(dist.cov)
    cov_names <- names(dist.cov)
    if(is.null(cov_names)) stop('dist.cov must be a named list.')

    if(is.null(loc.cov)){
      n_loc_cov <- 0
      loc.cov <- vector('list', n)
    } else {
      if(is(loc.cov, 'data.frame')){
        loc.cov <- list(loc.cov)
        n_loc_cov <- 1
      } else if(is(loc.cov, 'list')){
        n_loc_cov <- length(loc.cov)
      } else {
        stop('invlid input of the argument loc.cov.')
      }
    }

    if(!is.data.frame(mask)) {
      mask <- do.call('rbind', mask)
    }

    for(i in 1:n){
      cov_name <- cov_names[i]
      stopifnot(all(c('x', 'y') %in% colnames(dist.cov[[cov_name]])))

      tem <- dist.nearest(from = mask,
                         to = dist.cov[[cov_name]][, c('x', 'y')], col.name = cov_name)

      loc.cov[[n_loc_cov + i]] <- tem[!duplicated(tem[, c('x', 'y')]),]
    }

  }

  return(loc.cov)

}


convert.time.loc.cov.to.loc.cov <- function(time.loc.cov, loc.cov, session.cov){
  if(!is.null(time.loc.cov)){
    if(is.null(session.cov)) stop('please provide session.cov with column time')
    if(!all(c("session", "time") %in% colnames(session.cov))) stop('please provide information of session and time.')
    if(is(time.loc.cov, 'data.frame')) time.loc.cov <- list(time.loc.cov)
    if(!is(time.loc.cov, 'list')) stop('time.loc.cov should be a data frame or a list.')


    n <- length(time.loc.cov)

    if(is.null(loc.cov)){
      loc.cov <- vector('list', n)
      n_loc_cov <- 0
    } else {
      if(is(loc.cov, 'list')){
        n_loc_cov <- length(loc.cov)
      } else {
        loc.cov <- list(loc.cov)
        n_loc_cov <- 1
      }
    }

    session_time <- session.cov[,c('session', 'time'), drop = FALSE]
    for(i in 1:n){
      tem <- time.loc.cov[[i]]
      stopifnot(is(tem, 'data.frame'))
      stopifnot(all(c('time', 'x', 'y') %in% colnames(tem)))
      stopifnot(nrow(tem) > 3)
      u_time_from_session <- unique(session_time$time)
      u_time_from_loc_cov <- unique(tem$time)

      if(any(length(u_time_from_loc_cov) != length(u_time_from_session),
             !all(u_time_from_loc_cov %in% u_time_from_session))){
        stop('time provided in session.cov and time.loc.cov should match.')
      }

      tem <- merge(session_time, tem, by = 'time', all = T)

      tem <- tem[,-which(colnames(tem) == 'time'), drop = FALSE]

      loc.cov[[n_loc_cov + i]] <- tem[order(tem$session, tem$x, tem$y), ]
    }

    ## End of time.loc.cov.
  }

  return(loc.cov)
}

get.dataset.names <- function() {
  return(data(package="acre")$results[, "Item"])
}

separate.dist.loc.cov <- function(loc.cov, dist.cov.col.names) {
  ## Function that separates distance covariates from location covariates in a loc.cov list, which
  ## contains data frames of location and distance covariates.

  ## Check if loc.cov is a list.
  if (!is.list(loc.cov)) {
    stop("loc.cov must be a list of data frames.")
  }

  ## Check if dist.cov.col.names is a character vector.
  if (!is.character(dist.cov.col.names)) {
    stop("dist.cov.col.names must be a character vector of distance covariate column names.")
  }

  distance <- list()
  location <- list()

  for (i in seq_along(loc.cov)) {
    cov_df <- loc.cov[[i]]

    ## Check that the items of loc.cov are data frames.
    if (!is.data.frame(cov_df)) {
      stop("loc.cov must be a list of data frames.")
    }

    ## Make sure each data frame has x, y columns.
    if (!all(c("x", "y") %in% names(cov_df))) {
      stop("Each dataframe in loc.cov must contain 'x' and 'y' columns.")
    }

    ## Split the data frame into 2 new data frames, both with the same x, y
    ## columns, but one with only the distance covariates, and the other with only the location
    ## covariates.
    dist_cov_indices <- which(names(cov_df) %in% dist.cov.col.names)
    loc_cov_indices <- which(!names(cov_df) %in% c(dist.cov.col.names, "x", "y"))

    location_df <- cov_df[, c("x", "y", names(cov_df)[loc_cov_indices]), drop = FALSE]
    distance_df <- cov_df[, c("x", "y", names(cov_df)[dist_cov_indices]), drop = FALSE]

    ## Only add the data frame if it contains more than just the x, y columns.
    if (ncol(location_df) > 2) {
      location[[length(location) + 1]] <- location_df
    }

    if (ncol(distance_df) > 2) {
      distance[[length(distance) + 1]] <- distance_df
    }
  }

  return(list(distance = distance, location = location))
}

conf.int.column.can.be.non.scientific <- function(x) {
  all(abs(x) >= 1e-6 & abs(x) <= 1e6, na.rm = TRUE)
}

conf.int.matrix.format <- function(mat) {
  ## Format the matrix according to the desired behaviour:
  ##
  ## First two columns: if ALL values in that column can be printed as
  ## non-scientific, print as such for that column, else print in scientific.
  ##
  ## Second two columns: if ALL values ACROSS BOTH columns can be printed as
  ## non-scientific, print as such for both columns, else print in scientific.
  ##
  ## This is a little hacky; there is definitely a better way.

  formatted_mat <- mat
  ## Format the first two columns independently.
  for (col in 1:2) {
    if (conf.int.column.can.be.non.scientific(mat[, col])) {
      formatted_mat[, col] <- format(mat[, col], scientific = FALSE, trim=T)
    } else {
      formatted_mat[, col] <- format(mat[, col], scientific = TRUE, trim=T)
    }
  }

  ## Format the last two columns together.
  if (all(sapply(mat[, 3:4], conf.int.column.can.be.non.scientific))) {
    formatted_mat[, 3:4] <- format(mat[, 3:4], scientific = FALSE, trim=T)
  } else {
    formatted_mat[, 3:4] <- format(mat[, 3:4], scientific = TRUE, trim=T)
  }

  return(formatted_mat)
}

is.animal.model <- function(fit) {
  return(ifelse(is.null(fit$arg_input$captures$animal_ID), FALSE, TRUE))
}

faded.virdis <- function(n.levels, max.alpha = 1, min.alpha = 0) {
  ## Generate base viridis colours.
  base_colors <- viridisLite::viridis(n.levels)

  ## Generate a sequence of alpha scales.
  alpha_values <- seq(min.alpha, max.alpha, length.out = n.levels)

  ## Apply the alphas to our base colours.
  faded_colors <- unlist(
    Map(function(col, a) scales::alpha(col, a), base_colors, alpha_values)
  )

  return(faded_colors)
}
