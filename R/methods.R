#' Extract coefficients from the output of acre model
#'
#' @param object a fitted model from fit.acre().
#' @param types a character vector, accept any subset from ('all', 'fitted', 'linked', 'derived'), default is 'linked'. In details:
#'              linked - There is a link function attached to each parameter, either 'log', 'logit' or 'identical'. The 'linked' 
#'                       estimation is the estimation before back-transformation according to the link functions. For any extended
#'                       parameters, the estimations of the coefficients of any covariates will be always under link function of
#'                       identity. For example, if we have D ~ x1, the estimations of D_int and D_x1 are 1 and 0.5, if the value
#'                       of 'x1' is not provided in the argument of 'new.covariates', then no matter what the 'types' is, these
#'                       estimations will be shown all the time. If x1 = 3, then when 'linked', the function will return 1 + 0.5 * 3
#'                       as the estimation of log(D); when 'fitted', the function will return exp(1 + 0.5 * 3) as the estimation of D.
#'              fitted - The back transformed estimations from the 'linked' estimations and their link function.
#'              derived - The estimations of 'esa' for each session.
#'                       
#' @param pars a character vector containing any parameter names.
#' @param new.covariates a data frame containing the values of covariates of any extended parameter.
#' @param ... For S3 compatibility.
#'
#' @return a named numeric vector
#' @export
#'
coef.acre <- function(object, types = NULL, pars = NULL, new.covariates = NULL, ...){
  newdata <- new.covariates
  ## Deal with the default setting for 'types'.
  tem <- types.pars.sol(types, pars, newdata)
  ## tem = types.pars.sol("linked", pars, newdata)
  types <- tem$types
  pars <- tem$pars
  
  is.fitted <- 'fitted' %in% types
  is.linked <- 'linked' %in% types
  is.derived <- 'derived' %in% types
  
  ## Extract some key information from the object.
  param_values_og <- get.coef(object)
  name_og <- get.param.og(object)
  
  if (object$CL) {
    name_og <- setdiff(name_og, "D")
    pars <- setdiff(pars, "D")
  }
  
  name_extend <- get.par.extend.name(object)
  df_param <- get.data.param(object)
  esa <- esa(object)
  
  ## Check which parameters will be displayed.
  if(is.null(pars)){
    pars <- name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  if(is.null(name_extend) & !is.null(newdata)){
    warning('No parameter is extended, argument "newdata" will be ignored.')
    newdata <- NULL
  }
  
  output <- NULL
  
  ## Since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  ## consistent with the default order of all parameters.
  pars <- fulllist.par.generator()[fulllist.par.generator() %in% pars]
  extend_par_covariates <- get.par.extend.covariate(object)
  
  for(i in pars){
    par_info <- subset(df_param, df_param$par == i)
    values_link <- param_values_og[[i]]
    ## This indicates whether we calculate \beta \times X, i.e., for an extended parameter,
    ## whether we evaluate the link function.
    do_BX <- FALSE
    
    if(!i %in% name_extend){
      names(values_link) <- paste(i, "link", sep = "_")
    } else {
      names(values_link) <- gsub(" _ ", "\\.", names(values_link))
      names(values_link) <- paste(names(values_link), "link", sep = "_")
      if(all(extend_par_covariates[[i]] %in% colnames(newdata))){
        do_BX <- TRUE
      }
    }
    
    if(is.linked){
      if(!do_BX){
        output <- c(output, values_link)
      } else {
        gam <- get.gam(object, i)
        values_link_evaluated <- get.extended.par.value(gam, par_info$n_col_full, par_info$n_col_mask,
                                                       values_link, newdata)
        names(values_link_evaluated) <- rep(paste(i, "link", sep = "_"), length(values_link_evaluated))
        output <- c(output, values_link_evaluated)
      }

    }
    
    if(is.fitted){
      link <- par_info$link
      if(i %in% name_extend){
        if(!do_BX){
          ## If no new covariates are assigned and the parameter is extended, then treat all of the
          ## relevant betas as identity-linked parameters.
          values_fitted <- values_link
          names(values_fitted) <- gsub("_link", "", names(values_fitted))
          link <- 'identity'
        } else {
          gam <- get.gam(object, i)
          values_fitted <- get.extended.par.value(gam, par_info$n_col_full,
                                                 par_info$n_col_mask, values_link, newdata)

          names(values_fitted) <- rep(i, length(values_fitted))

        }
      } else {
        values_fitted <- values_link
        names(values_fitted) <- i
      }

      values_fitted <- unlink.fun(link = link, value = values_fitted)
      output <- c(output, values_fitted)
    }
  }
  
  if(is.derived){
    value_esa <- esa$value
    names(value_esa) <- paste('esa', seq(length(value_esa)), sep = '.')
    output <- c(output, value_esa)
  }
  
  class(output) <- 'coef.acre'
  
  return(output)
}


#' Title
#'
#' @param object an object generated by the bootstrap function "boot.acre()".
#' @param types a character vector, the same as "coef.acre()".
#' @param pars a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param correct.bias logical. if TRUE, apply a bias correction method to the bootstrap results; otherwise,
#'                     return the original estimation directly. The bias correct method is
#'                     hat(coef) - (mean(boot_results) - hat(coef)).
#' @param ... For S3 compatibility.
#' @export
coef.acreboot <- function(object, types = NULL, pars = NULL, new.covariates = NULL, correct.bias = FALSE, ...){
  
  if(!correct.bias){
    output <- coef.acre(object = object, types = types, pars = pars, new.covariates = new.covariates)
    return(output)
  }

  
  ## Some preparation:
  ## deal with the default setting for 'types'.
  tem <- types.pars.sol(types, pars, new.covariates)
  types <- tem$types
  pars <- tem$pars
  name_og <- get.param.og(object)

  ## Check which parameters will be displayed.
  if(is.null(pars)){
    pars <- name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  ## Since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  ## consistent with the default order of all parameters.
  pars <- fulllist.par.generator()[fulllist.par.generator() %in% pars]
  
  ## Get the basic information needed for this method.
  name_extend <- get.par.extend.name(object)
  
  if(is.null(name_extend) & !is.null(new.covariates)){
    warning('No parameter is extended, argument "new.covariates" will be ignored.')
    new.covariates <- NULL
  }
  
  

  df_param <- get.data.param(object)
  res <- get.boot.res(object, pars)
  coefs <- coef.acre(object, types = 'linked', pars)
  bias <- get.bias(res, coefs)
  coefs_linked_bias_corrected <- coefs - bias
  
  linked_name <- names(coefs)
  original_name <- ori.name(linked_name)
  extend_par_covariates <- get.par.extend.covariate(object)
  
  output <- vector('list', length(types))
  names(output) <- types
  
  for(i in types){
    
    if(i == 'linked'){
      if(is.null(new.covariates)){
        output[[i]] <- coefs_linked_bias_corrected
      } else {
        ## If new.covariates is provided, then we need to go through the parameters one by one.
        for(j in pars){
          values_link <- coefs_linked_bias_corrected[original_name == j]
          par_info <- subset(df_param, df_param$par == j)
          
          if(j %in% name_extend){
            if(!all(extend_par_covariates[[j]] %in% colnames(new.covariates))){
              values_link_evaluated <- values_link
            } else {
              ## if(nrow(new.covariates) != 1) stop('Argument "new.covariates" can only accept 1 row.')
              gam <- get.gam(object, j)
              values_link_evaluated <- get.extended.par.value(gam, par_info$n_col_full, par_info$n_col_mask,
                                                     values_link, new.covariates)
              
              names(values_link_evaluated) <- rep(paste(j, 'link', sep = '_'), length(values_link_evaluated))
            }
          } else {
            values_link_evaluated <- values_link
          }
          
          output[[i]] <- c(output[[i]], values_link_evaluated)
        }
        
      }
      
    } else if (i == 'derived'){
      coefs_esa <- coef.acre(object, types = 'derived')
      res_esa <- get.boot.res.esa(object)
      bias_esa <- get.bias(res_esa, coefs_esa)
      output[[i]] <- coefs_esa - bias_esa
      
    } else {
      ## When types is 'fitted', use 'coefs_linked_bias_corrected' and new.covariates to
      ## calculate the estimates.

      for(j in pars){
        values_link <- coefs_linked_bias_corrected[original_name == j]
        
        par_info <- subset(df_param, df_param$par == j)
        link <- par_info$link
        
        if(j %in% name_extend){
          if(!all(extend_par_covariates[[j]] %in% colnames(new.covariates))){
            ## If no new covariates are assigned and the parameter is extended, then treat all of the
            ## relevant betas as identity-linked parameters.
            values_fitted <- values_link
            names(values_fitted) <- gsub("_link", "", names(values_fitted))
            link <- 'identity'
          } else {
            ## if(nrow(new.covariates) != 1) stop('Argument "new.covariates" can only accept 1 row.')
            
            gam <- get.gam(object, j)
            values_fitted <- get.extended.par.value(gam, par_info$n_col_full, par_info$n_col_mask,
                                                   values_link, new.covariates)
            
            names(values_fitted) <- rep(j, length(values_fitted))
          }
        } else {
          values_fitted <- values_link
          names(values_fitted) <- j
        }
        
        values_fitted <- unlink.fun(link = link, value = values_fitted)

        output[[i]] <- c(output[[i]], values_fitted)
      }
      ## End of if(i == 'fitted').
    }
    
    ## End of for i in types.
  }

  names(output) <- NULL
  output <- do.call('c', output)

  
  class(output) <- "coef.acre"
  return(output)
  
}

#' @export
print.coef.acre <- function(x, ...){
  output <- as.matrix(x)
  colnames(output) <- "Est"
  print(output)
}

#' Extract variance covariance matrix of the estimated parameters from acre models
#'
#' @param object a fitted model from fit.acre().
#' @param types a character vector, the same as "coef.acre()".
#' @param pars a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param show.fixed.par a logical value. To control whether to include the fixed parameters in the covariance matrix.
#'                       It is TRUE by default.
#' @param ... For S3 compatibility.
#'
#' @return a list with matrices as its elements if multiple 'types'. a matrix if only one 'types'.
#' @export
vcov.acre <- function(object, types = NULL, pars = NULL, new.covariates = NULL, show.fixed.par = TRUE, ...){
  
  ## Deal with the default setting for 'types'.
  tem <- types.pars.sol(types, pars, new.covariates)
  types <- tem$types
  pars <- tem$pars
  
  ## Below, 'og' means the original output from the object.
  param_values_og <- get.coef(object)
  cov_og <- object$vcov
  name_dim_og <- colnames(cov_og)
  name_dim <- gsub(" _ ", "\\.", name_dim_og)
  which.derived <- which(substr(name_dim, 1, 3) == 'esa')
  name_dim <- name_dim[-which.derived]
  name_dim_og <- name_dim_og[-which.derived]
  cov_derived <- cov_og[which.derived, which.derived, drop = FALSE]
  cov_linked <- cov_og[-which.derived, -which.derived, drop = FALSE]
  dimnames(cov_linked) <- list(name_dim, name_dim)
  name_extend <- get.par.extend.name(object)
  extend_par_covariates <- get.par.extend.covariate(object)
  
  
  if(is.null(name_extend) & !is.null(new.covariates)){
    warning('No parameter is extended, argument "new.covariates" will be ignored.')
    new.covariates <- NULL
  }
  name_extned_covariate_provided <- ext.par.in.new.df(name_extend, new.covariates, extend_par_covariates)
  
  ## The original name is needed: in "name_og", "D.(Intercept)", "D.x1" and "D.xxx" are all regarded as 'D' only,
  ## and the same applies to the other extendable parameters.
  
  ## Here 'og' means the original parameter names.
  name_og <- ori.name(name_dim)
  
  ## Extract the link functions for all elements in the linked covariance matrix;
  ## if new.covariates is provided, this "link.funs" will be useless, as the covariance matrices of the fitted and
  ## linked parameters will have different dimensions.
  df_param <- get.data.param(object)
  link_funs <- character(length(name_og))
  ## Extract the linked estimated values. The reason for not using do.call('c', param_values_og)
  ## is that the fixed parameters could mess this up.
  param_values <- numeric(length(name_og))
  
  for(i in 1:length(name_dim)){
    param_values[i] <- param_values_og[[name_og[i]]][name_dim_og[i]]
    if(!name_og[i] %in% name_extend){
      name_dim[i] <- name_og[i]
      link_funs[i] <- df_param[which(df_param$par == name_og[i]), 'link']
    } else {
      link_funs[i] <- 'identity'
    }
  }
  
  ## Check which parameters will be displayed.
  if(is.null(pars)){
    pars <- unique(name_og)
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters, which is not fixed, in this model.")
  }
  
  
  ## Since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  ## consistent with the default order of all parameters.
  pars <- fulllist.par.generator()[fulllist.par.generator() %in% pars]
  ## The vector name_og is naturally in the right order because it is extracted from the TMB variance-covariance matrix output.
  
  ## Obtain the indices for the assigned "pars".
  index_par <- which(name_og %in% pars)
  name_og <- name_og[index_par]
  name_dim <- name_dim[index_par]
  cov_linked <- cov_linked[index_par, index_par, drop = FALSE]
  link_funs <- link_funs[index_par]
  param_values <- param_values[index_par]
  fixed_par <- get.fixed.par.name(object)
  
  ## Suppress "D" if the conditional likelihood was fitted.
  if (object$CL) {
    fixed_par <- setdiff(fixed_par, "D")
  }
  
  output <- vector('list', length(types))
  names(output) <- types
  for(type in types){
    if(type == "linked"){
      if(is.null(new.covariates)){
        output[[type]] <- cov_linked
        tem <- paste(name_dim, 'link', sep = '_')
        dimnames(output[[type]]) <- list(tem, tem)
      } else {
        gam.output <- get.gam(object)
        output[[type]] <- delta.method.acre(cov_linked, param_values, new.covariates = new.covariates,
                                               name.og = name_og, name.extend = name_extend,
                                               par.ext.cov.provided = name_extned_covariate_provided, df.param = df_param,
                                               gam.output = gam.output, back.trans = FALSE)
        
        new_name <- character(0)
        for(p in pars){
          if(p %in% name_extend){
            if(p %in% name_extned_covariate_provided){
              tem <- rep(paste0(p, "_link"), nrow(new.covariates))
            } else {
              tem <- paste(name_dim[which(name_og == p)], "link", sep = "_")
            }
          } else {
            tem <- paste0(p, "_link") 
          }
          new_name <- c(new_name, tem)
        }
        
        dimnames(output[[type]]) <- list(new_name, new_name)
        
      }
      
      
      if(show.fixed.par){
        output[[type]] <- fixed.par.vcov.add(output[[type]], fixed_par, type)
      }

    } else if(type == 'derived'){
      output[[type]] <- cov_derived
    } else if(type == 'fitted'){
      ## If 'new.covariates' is not provided, just keep the extended covariates as parameters with an identity link function.
      if(is.null(new.covariates)){
        output[[type]] <- delta.method.acre(cov_linked, param_values, link.funs = link_funs, back.trans = TRUE)
        dimnames(output[[type]]) <- list(name_dim, name_dim)
      } else {
        gam.output <- get.gam(object)
        output[[type]] <- delta.method.acre(cov_linked, param_values, new.covariates = new.covariates,
                                               name.og = name_og, name.extend = name_extend,
                                               par.ext.cov.provided = name_extned_covariate_provided, df.param = df_param,
                                               gam.output = gam.output, back.trans = TRUE)
        
        new_name <- character(0)
        for(p in pars){
          if(p %in% name_extend){
            if(p %in% name_extned_covariate_provided){
              tem <- rep(p, nrow(new.covariates))
            } else {
              tem <- name_dim[which(name_og == p)]
            }
          } else {
            tem <- p 
          }
          new_name <- c(new_name, tem)
        }
        
        dimnames(output[[type]]) <- list(new_name, new_name)
      }
      
      if(show.fixed.par){
        output[[type]] <- fixed.par.vcov.add(output[[type]], fixed_par, type)
      }
      
    }
  }
  
  ## If only one type is assigned, change the output from a list to a matrix.
  if(length(types) == 1) output <- output[[types]]
  
  return(output)
  
}


#' Title
#'
#' @param object an object generated by the bootstrap function "boot.acre()".
#' @param types a character vector, the same as "coef.acre()".
#' @param pars a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param from.boot a logical value, TRUE by default, to control whether the covariance matrix is calculated
#'                  from the bootstrap results. If FALSE, the method "vcov.acre()" will be called.
#' @param show.fixed.par a logical value, the same as "vcov.acre()".
#' @param ... For S3 compatibility.
#' @export
vcov.acreboot <- function(object, types = NULL, pars = NULL, new.covariates = NULL, from.boot = TRUE, show.fixed.par = TRUE, ...){
  
  if(!from.boot){
    output <- vcov.acre(object = object, types = types, pars = pars, new.covariates = new.covariates, show.fixed.par = show.fixed.par)
    return(output)
  }
  
  ## Deal with the default setting for 'types'.
  tem <- types.pars.sol(types, pars, new.covariates)
  types <- tem$types
  pars <- tem$pars
  all_par_name_og <- get.param.og(object)
  ## Check which parameters will be displayed.
  if(is.null(pars)){
    pars <- all_par_name_og
  } else {
    if(any(!pars %in% all_par_name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  ## Since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  ## consistent with the default order of all parameters.
  pars <- fulllist.par.generator()[fulllist.par.generator() %in% pars]
  if(any(c('linked', 'fitted') %in% types)){
    res <- get.boot.res(object, pars)
    fixed_par <- get.fixed.par.name(object)
    pars <- pars[!pars %in% fixed_par]
  }
  output <- vector('list', length(types))
  names(output) <- types
  for(i in types){
    if(i == 'linked'){
      res_linked <- res.transform(res, new.covariates, pars, object, back.trans = FALSE)
      output[[i]] <- var.from.res(res_linked, fixed_par)
      
      if(show.fixed.par){
        output[[i]] <- fixed.par.vcov.add(output[[i]], fixed_par, i)
      }
      
    } else if(i == 'derived'){
      res_esa <- get.boot.res.esa(object)
      output[[i]] <- var.from.res(res_esa)
    } else {
      res_fitted <- res.transform(res, new.covariates, pars, object, back.trans = TRUE)
      output[[i]] <- var.from.res(res_fitted, fixed_par)
      
      if(show.fixed.par){
        output[[i]] <- fixed.par.vcov.add(output[[i]], fixed_par, i)
      }
      
    }
  }
  
  ## When there is only one type, output it directly without a list structure.
  if(length(types) == 1) output <- output[[types]]
 
  return(output)

    
}


#' Title
#'
#' @param object 
#' @param ... 
#' @export
stdEr <- function(object, ...){
  UseMethod("stdEr")
}




##################################################################################################################

#' Extract standard errors of the estimated parameters from acre models
#'
#' @param object a fitted model from "fit.acre()".
#' @param types a character vector, the same as "coef.acre()".
#' @param pars a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param show.fixed.par a logical value, the same as "vcov.acre()".
#' @param ... For S3 compatibility.
#'
#' @return a named numeric vector
#' @export
stdEr.acre <- function(object, types = NULL, pars = NULL, new.covariates = NULL, show.fixed.par = TRUE, ...){
  
  output_vcov <- vcov.acre(object = object, types = types, pars = pars, new.covariates = new.covariates,
                              show.fixed.par = show.fixed.par)
  if(is(output_vcov, 'list')){
    output <- NULL
    for(i in names(output_vcov)){
      output <- c(output, sqrt(diag(output_vcov[[i]])))
    }
  } else {
    output <- sqrt(diag(output_vcov))
  }
  class(output) <- 'std.acre'
  return(output)
}


#' Title
#'
#' @param object an object generated by the bootstrap function "boot.acre()".
#' @param types a character vector, the same as "coef.acre()".
#' @param pars a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param from.boot a logical value, the same as "vcov.acreboot()".
#' @param show.fixed.par a logical value, the same as "vcov.acre()".
#' @param ... For S3 compatibility.
#' @export
stdEr.acreboot <- function(object, types = NULL, pars = NULL, new.covariates = NULL, from.boot = TRUE, show.fixed.par = TRUE, ...){

  output_vcov <- vcov.acreboot(object = object, types = types, pars = pars, new.covariates = new.covariates,
                               from.boot = from.boot, show.fixed.par = show.fixed.par)
  
  if(is(output_vcov, 'list')){
    output <- NULL
    for(i in names(output_vcov)){
      output <- c(output, sqrt(diag(output_vcov[[i]])))
    }
  } else {
    output <- sqrt(diag(output_vcov))
  }
  class(output) <- 'std.acre'
  return(output)
  
  
}

#' @export
print.std.acre <- function(x, ...){
  output <- as.matrix(x)
  colnames(output) <- "Std"
  print(output)
}

######################################################################################################################

#' Extract confidence interval for acre.tmb models
#'
#' @param object a fitted model from "fit.acre()".
#' @param parm a character vector, the same as "coef.acre()".
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param types a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param ... For S3 compatibility.
#'
#' @return a matrix
#' @export
confint.acre <- function(object, parm = NULL, level = 0.95, types = NULL, new.covariates = NULL, ...){
  ## Deal with the default setting for 'types'.
  tem <- types.pars.sol(types, parm, new.covariates)
  types <- tem$types
  pars <- tem$pars

  
  stopifnot(all(level < 1 & level > 0))
  p_upper <- 0.5 + 0.5 * level
  p_lower <- 0.5 - 0.5 * level
  col_name <- paste(c(p_lower, p_upper) * 100, "%", sep = " ")
  q_upper <- qnorm(p_upper)
  q_lower <- qnorm(p_lower)

  output <- gaussian.confint.cal(object = object, types = types, pars = pars,
                                new.covariates = new.covariates, q.lower = q_lower,
                                q.upper = q_upper)

  df_param <- get.data.param(object)
  if('fitted' %in% types){
    tem <- output[['fitted']]
    name_og <- get.param.og(object)
    for(i in 1:nrow(tem)){
      p <- tem[i, 'par']
      if(p %in% name_og){
        link <- df_param[which(df_param$par == p), 'link']
      } else {
        link <- 'identity'
      }
      
      tem[i, c('lower', 'upper')] <- unlink.fun(link, tem[i, c('lower', 'upper')])
    }
    output[['fitted']] <- tem
  }
  
  for(i in names(output)){
    par_name <- output[[i]]$par
    output[[i]] <- as.matrix(output[[i]][, c('lower', 'upper'), drop = FALSE])
    rownames(output[[i]]) <- par_name
  }
  
  output <- do.call('rbind', output)
  colnames(output) <- col_name
  return(output)
}




#' Title
#'
#' @param object an object generated by the bootstrap function "boot.acre()".
#' @param parm a character vector, the same as "coef.acre()".
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param types a character vector, the same as "coef.acre()".
#' @param new.covariates a data frame, the same as "coef.acre()".
#' @param correct.bias a logical value indicates whether to apply a bias correction method to
#'                     the confidence interval. Default is FALSE, and the naive percentiles from
#'                     the bootstrap results will used as confidence interval; if TRUE, the 
#'                     bootstrap results will be corrected by beta_boot = 2 * beta_est - beta_boot,
#'                     where beta_boot is the bootstrap results values and beta_est is the
#'                     estimations from "coef.acre()", and the percentiles from these corrected
#'                     bootstrap results are used as the confidence interval.
#' @param from.boot a logical value indicates whether to use the bootstrap results to construct 
#'                  the confidence matrix. Default is TRUE; if FALSE, "confint.acre()" will
#'                  be called.
#' @param ... For S3 compatibility.
#' @export
confint.acreboot <- function(object, parm = NULL, level = 0.95, types = NULL, new.covariates = NULL,
                             correct.bias = FALSE, from.boot = TRUE, ...){
  if(!from.boot){
    if(correct.bias) message(paste0("The argument 'correct.bias' will be ignored as the confidence",
      " interval is not from bootstrap results."))
    
    output <- confint.acre(object = object, types = types, level = level, pars = parm,
                              new.covariates = new.covariates)
    return(output)
  }
  
  stopifnot(all(level < 1 & level > 0))
  name_og <- get.param.og(object)
  tem <- types.pars.sol(types, parm, new.covariates)
  types <- tem$types
  pars <- tem$pars
  
  
  ## Check which parameters will be displayed.
  if(is.null(pars)){
    pars <- name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  ## Make sure pars are in the right order.
  pars <- fulllist.par.generator()[fulllist.par.generator() %in% pars]

  
  ## Obtain the bootstrap result.
  res <- get.boot.res(object, pars)
  coef_link <- coef.acre(object, types = 'linked', pars)
  
  res <- res.mod.for.CI(res, coef_link, correct.bias)

  
  ## Prepare for output.
  output <- vector('list', length(types))
  ## Deal with different output types one by one.
  for(i in types){
    if(i == 'linked'){
      res_linked <- res.transform(res, new.covariates, pars, object, back.trans = FALSE)
      output[[i]] <- boot.res.to.CI(res_linked, level)
    } else if(i == 'derived'){
      res_esa <- get.boot.res.esa(object)
      coef_esa <- coef.acre(object, types = 'derived')
      res_esa <- res.mod.for.CI(res_esa, coef_esa, correct.bias)
      
      output[[i]] <- boot.res.to.CI(res_esa, level)
    } else {
      res_fitted <- res.transform(res, new.covariates, pars, object, back.trans = TRUE)
      output[[i]] <- boot.res.to.CI(res_fitted, level)
      
    }
  }
  
  output <- do.call('rbind', output)

  return(output)
}



######################################################################################################################

#' Calculate AIC for acre Models
#'
#' @param ... One or more acre model objects.
#' @param k Penalty parameter (default is 2 for AIC).
#'
#' @return A named vector of AIC values or NA if conditions are not met.
#' @export
AIC.acre <- function(..., k = 2) {
  ## Capture all model objects.
  models <- list(...)
  
  ## Check if all objects are of class 'acre'.
  if (!all(sapply(models, inherits, what = "acre"))) {
    stop("All inputs must be of class 'acre'.")
  }
  
  ## Initialise result storage.
  aic_values <- numeric(length(models))
  n_par <- numeric(length(models))
  names(aic_values) <- match.call(expand.dots = FALSE)$`...`
  
  ## Loop over models and calculate AIC.
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    ## A call model fitted with known cue rates.
    fit.freqs <- !is.null(model$args$cue.rates) && !"animal_ID" %in% colnames(model$args$capt)
    if (fit.freqs) {
      message(sprintf("NOTE: AIC for model '%s' relies on independence between locations of calls from the same animal, which may not be appropriate.", names(aic_values)[i]))
      aic_values[i] <- NA
    } else {
      n_fixed <- length(get.fixed.par.name(model))
      aic_values[i] <- k * (length(coef(model)) - n_fixed) - 2 * model$loglik
      n_par[i] <- length(coef(model))
    }
  }
  
  return(data.frame(npar = n_par, AIC = aic_values))
}



#' Title
#'
#' @param object a fitted model from "fit.acre()".
#' @param type a character vector, could be either "response" or "link". The "response" correspondence to
#'             the type of "fitted" in the function "coef.acre()" and the "link" correspondence to
#'             the type of "linked" in the function "coef.acre()".
#' @param newdata a data frame, the same as "coef.acre()". If there is any extended parameters, and
#'                the corresponding covariates is not provided here, the estimation for such parameters
#'                will not be shown.
#' @param se.fit a logical value indicates whether to show standard error.
#' @param confidence a logical value indicates whether to show confidence interval.
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param realnames a character vector containing any parameter names.
#' @param ... For S3 compatibility.
#' @export
predict.acre <- function(object, type = 'response', newdata = NULL, se.fit = TRUE, confidence = TRUE,
                            level = 0.95, realnames = NULL, ...){
  fit <- object
  stopifnot(all(type %in% c('link', 'response'), length(type) == 1))
  
  if(type == 'response'){
    type <- 'fitted'
  } else if(type == 'link'){
    type <- 'linked'
  } else {
    stop('invalid "type", it could be either "response" or "link".')
  }
  
  name_og <- get.param.og(fit)
  
  ## Check which parameters will be displayed.
  if(is.null(realnames)){
    realnames <- name_og
  } else {
    if(any(!realnames %in% name_og)) stop("Argument 'realnames' only accept parameters' name in this model.")
  }
  
  if (fit$CL) {
    realnames <- setdiff(realnames, "D")
  }
  
  ## Make sure pars are in the right order.
  realnames <- fulllist.par.generator()[fulllist.par.generator() %in% realnames]
  
  par_fixed_name <- get.fixed.par.name(fit)
  par_ext_name <- get.par.extend.name(fit)
  par_ext_vars <- get.par.extend.covariate(fit)
  newdata_vars <- colnames(newdata)
  
  output <- vector('list', length(realnames))
  names(output) <- realnames
  
  for(i in realnames){
    if(i %in% par_ext_name){
      if(!all(par_ext_vars[[i]] %in% newdata_vars)){
        n_row <- 0
        output[[i]] <- paste0('Estimates not provided, because ', i, ' depends on ', paste(par_ext_vars[[i]], collapse = ', '),
                             ', which are not all provided in the argument "newdata".')
      } else {
        ## If the condition above is FALSE, then newdata must not be NULL.
        n_row <- nrow(newdata)
      }
      
    } else {
      n_row <- 1
    }
    
    if(n_row != 0) {
      output[[i]] <- data.frame(Estimate = numeric(n_row))

      output[[i]]$Estimate <- as.vector(coef(fit, types = type, new.covariates = newdata, pars = i))
      
      if(se.fit){
        if(i %in% par_fixed_name){
          output[[i]]$StdError <- 0
        } else {
          output[[i]]$StdError <- as.vector(stdEr(fit, types = type, new.covariates = newdata, pars = i, show.fixed.par = FALSE))
        }
        
      } 
      if(confidence){
        if(i %in% par_fixed_name){
          output[[i]]$Lower <- output[[i]]$Estimate
          output[[i]]$Upper <- output[[i]]$Estimate
        } else {
          output[[i]]$Lower <- numeric(n_row)
          output[[i]]$Upper <- numeric(n_row)
          output[[i]][, c('Lower', 'Upper')] <- confint(fit, parm = i, types = type, level = level, new.covariates = newdata)
        }

      }
    }

  }
  
  class(output) <- 'predict.acre'
  return(output)
  
}



#' Title
#'
#' @param object an object generated by the bootstrap function "boot.acre()".
#' @param type a character vector, the same as "predict.acre()".
#' @param newdata a data frame, the same as "predict.acre()".
#' @param se.fit a logical value indicates whether to show standard error.
#' @param confidence a logical value indicates whether to show confidence interval.
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param realnames a character vector, the same as "predict.acre()".
#' @param correct.bias a logical value indicates whether to apply the bias correction method to the
#'                     point estimations and the confidence intervals. The method is described in the
#'                     functions "coef.acreboot()" and "confint.acreboot()".
#' @param from.boot a logical value indicates whether to extract standard error and confidence interval
#'                  from the bootstrap results.
#' @param ... For S3 compatibility.
#' @export
predict.acreboot <- function(object, type = 'response', newdata = NULL, se.fit = TRUE, confidence = TRUE,
                             level = 0.95, realnames = NULL, correct.bias = FALSE, from.boot = TRUE, ...){
  fit <- object
  if(!from.boot){
    if(correct.bias) message(paste0("The argument 'correct.bias' will be ignored as the confidence",
                                    " interval is not from bootstrap results."))
    correct.bias <- FALSE
  }
  
  
  stopifnot(all(type %in% c('link', 'response'), length(type) == 1))
  
  if(type == 'response'){
    type <- 'fitted'
  } else {
    type <- 'linked'
  }
  
  name_og <- get.param.og(fit)
  
  ## Check which parameters will be displayed.
  if(is.null(realnames)){
    realnames <- name_og
  } else {
    if(any(!realnames %in% name_og)) stop("Argument 'realnames' only accept parameters' name in this model.")
  }
  
  ## Make sure pars are in the right order.
  realnames <- fulllist.par.generator()[fulllist.par.generator() %in% realnames]
  
  par_fixed_name <- get.fixed.par.name(fit)
  par_ext_name <- get.par.extend.name(fit)
  par_ext_vars <- get.par.extend.covariate(fit)
  newdata_vars <- colnames(newdata)
  
  output <- vector('list', length(realnames))
  names(output) <- realnames
  
  for(i in realnames){
    if(i %in% par_ext_name){
      if(!all(par_ext_vars[[i]] %in% newdata_vars)){
        n_row <- 0
        output[[i]] <- paste0('Estimates not provided, because ', i, ' depends on ', paste(par_ext_vars[[i]], collapse = ', '),
                             ', which are not all provided in the argument "newdata".')
      } else {
        ## If the condition above is FALSE, then newdata must not be NULL.
        n_row <- nrow(newdata)
      }
      
    } else {
      n_row <- 1
    }
    
    if(n_row != 0){
      output[[i]] <- data.frame(Estimate = numeric(n_row))

      output[[i]]$Estimate <- as.vector(coef(fit, types = type, new.covariates = newdata, pars = i, correct.bias = correct.bias))
      
      if(se.fit){
        
        if(i %in% par_fixed_name){
          output[[i]]$StdError <- 0
        } else {
          output[[i]]$StdError <- as.vector(stdEr(fit, types = type, new.covariates = newdata, pars = i, from.boot = from.boot, show.fixed.par = FALSE))
        }
        
      }

      if(confidence){
        if(i %in% par_fixed_name){
          output[[i]]$Lower <- output[[i]]$Estimate
          output[[i]]$Upper <- output[[i]]$Estimate
        } else {
          output[[i]]$Lower <- numeric(n_row)
          output[[i]]$Upper <- numeric(n_row)
          output[[i]][, c('Lower', 'Upper')] <- confint(fit, parm = i, types = type, new.covariates = newdata, correct.bias = correct.bias,
                                                       from.boot = from.boot)
        }
      }
    }
    
  }
  
  class(output) <- 'predict.acre'
  return(output)
  
}

#' @export
print.predict.acre <- function(x, ...){
  for(i in names(x)){
    cat(paste0(i, ": \n"))
    if(is(x[[i]], 'data.frame')){
      print(x[[i]], row.names = FALSE)
      cat("\n")
    } else {
      cat(paste0(x[[i]], "\n"))
      cat("\n")
    }
  }
}

#' Summarise `acre` Model Fits
#'
#' @param object an object generated from the model fitting function "fit.acre()" or
#'               the bootstrap process "boot.acre()".
#' @param derived.print a logical value indicates whether to show the estimations of "esa".
#' @param ... For S3 compatibility.
#' @export
summary.acre <- function(object, derived.print = FALSE, ...){
  coefs <- coef(object, types = 'fitted')
  derived <- coef(object, types = 'derived')
  coefs_se <- stdEr(object, types = 'fitted', show.fixed.par = FALSE)
  derived_se <- stdEr(object, types = 'derived', show.fixed.par = FALSE)
  CI <- confint(object, types = 'fitted')
  CI_derived <- confint(object, types = 'derived')
  is_boot <- is(object, 'acreboot')
  ## Conditional likelihood adjustment to the output (remove "D").
  CL <- object$CL
  if (object$CL) {
    coefs <- coefs[!grepl("^D", names(coefs))]
  }
  
  infotypes <- get.infotypes(object)
  detfn <- get.detfn(object)
  n.sessions <- object$n.sessions
  
  pars_ext <- get.par.extend.name(object)
  
  if(!is.null(pars_ext)){
    df_parm <- get.data.param(object)
    df_parm <- df_parm[df_parm$par %in% pars_ext, , drop = FALSE]
    pars_ext_links <- df_parm$link
    names(pars_ext_links) <- df_parm$par
  } else {
    pars_ext_links <- NULL
  }
  
  if('ss' %in% infotypes){
    ss_link <- get.ss.link(object)
    ss_cutoff <- get.ss.opts(object)$cutoff
    ss_opts <- list(cutoff = ss_cutoff, link = ss_link)
  } else {
    ss_opts <- NULL
  }
  
  output <- list(coefs = coefs, derived = derived,
                coefs_se = coefs_se, derived_se = derived_se,
                CI = CI, CI_derived = CI_derived,
                is_boot = is_boot, infotypes = infotypes, detfn = detfn,
                n.sessions = n.sessions, pars_ext_links = pars_ext_links,
                ss_opts = ss_opts, derived_print = derived.print, CL = CL)
  class(output) <- c("summary.acre", class(output))
  return(output)
}

#' @export
print.summary.acre <- function(x, ...){
  
  derived_print <- x$derived_print
  
  n.coefs <- length(x$coefs)
  n.derived <- length(x$derived)
  mat <- matrix(0, nrow = n.coefs + n.derived + 1, ncol = 4)
  mat[1:n.coefs, 1] <- as.vector(x$coefs)
  for(i in 1:n.coefs){
    name_coef <- names(x$coefs)[i]
    index_coefs_se <- which(names(x$coefs_se) == name_coef)
    if(length(index_coefs_se) != 0){
      mat[i, 2] <- as.vector(x$coefs_se[index_coefs_se])
    } else {
      mat[i, 2] <- 0
    }
  }
  
  mat[1:n.coefs, 3:4] <- x$CI
  mat[n.coefs + 1, ] <- NA
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), 1:2] <- cbind(x$derived, x$derived_se)
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), 3:4] <- x$CI_derived
  rownames(mat) <- c(names(x$coefs), "---", names(x$derived))
  colnames(mat) <- c("Estimate", "Std. Error", colnames(x$CI))
  
  if(!derived_print){
    mat <- mat[1:n.coefs,,drop = FALSE]
  }
  
  detfn <- c(hn = "Halfnormal", hhn = "Hazard halfnormal", hr = "Hazard rate", th = "Threshold",
             lth = "Log-link threshold", ss = "Signal strength")[x$detfn]
  infotypes <- c(bearing = "Bearings", dist = "Distances", ss = "Signal strengths",
                 toa = "Times of arrival", mrds = "Exact locations")[x$infotypes]
  
  pars_ext_links <- x$pars_ext_links
  ss_opts <- x$ss_opts
  cutoff <- ss_opts$cutoff
  ss_link <- ss_opts$link
  CI_method <- ifelse(x$is_boot, "Percentile", "Wald")
  
  cat("Detection function:", detfn, "\n")
  cat("Number of sessions:", x$n.sessions, "\n")
  cat("Information types: ")

  if(length(infotypes) != 0){
    cat(infotypes, sep = ", ")
  } else {
    cat("NULL")
  }
  
  cat("\n")
  cat("Confidence interval method:", CI_method, "\n")
  
  if (x$CL) {
    cat("Likelihood mode:", "Conditional", "\n")
  }
  
  cat("\n", "\n", "Parameters:", "\n")

  printCoefmat(mat, cs.ind = c(1, 3, 4), tst.ind = NULL)
  
  if(!is.null(pars_ext_links)){
    cat("\n", "\n", "Extended parameters link functions:", "\n")
    for(i in names(pars_ext_links)){
      cat(i, ":", pars_ext_links[i], "\n")
    }
  }
  
  if(!is.null(ss_opts)){
    cat("\n", "Signal strength related information:", "\n")
    cat("cutoff: ", cutoff, "\n")
    cat("link function:", ss_link, "\n")
  }
}


#' Extract (negative) log-likelihood for a fitted acre object
#'
#' @param object a fitted model from "fit.acre()".
#' @param ... For S3 compatibility.
#'
#' @return a numeric value, the negative log-likelihood of the model
#' @export
logLik.acre <- function(object, ...) {
  return(object$loglik)
}

#' Estimated effective sampling areas
#'
#' Extracts an estimate and standard error of each session's effective
#' sampling area.
#' 
#' @param fit A fitted model from [fit.acre()]
#' @export
esa <- function(fit){
  fit$esa
}

#' Estimated density at each mask point
#'
#' Extracts the estimated animal (or call) density at each mask point, for each session.
#'
#' @param fit A fitted model from [fit.acre()].
#' @param session The session to extract. If `NULL`, the default, all sessions are returned.
#'
#' @return For one session, a vector with the estimated density at each mask point; otherwise, a
#'         list of such vectors, one for each session.
#' @export
mask.density <- function(fit, session = NULL){
  if(is.null(session)){
    fit$D.mask
  } else {
    fit$D.mask[[session]]
  }
}

#' @export
print.acre <- function(x, ...){
  cat("Acoustic spatial capture-recapture model\n")
  if(!is.null(x$call)){
    cat("\nCall:\n")
    print(x$call)
  }
  est <- coef(x, types = "fitted")
  se <- stdEr(x, types = "fitted")
  cat("\nEstimates:\n")
  print(cbind(Estimate = unclass(est), "Std. Error" = unclass(se)[names(est)]), ...)
  cat("\nLog-likelihood:", format(x$loglik), "\n")
  invisible(x)
}

















