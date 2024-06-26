#' Bootstrapping a fitted `acre` model
#'
#' @param fit a model object generated by `fit.acre`.
#' @param N a numerical value indicates the number of iterations of the bootstrap process.
#' @param n.cores a numerical value indicates the number of cores to be used. Currently this function doesn't support
#'                parallel computing, so it is fixed to be 1.
#' @param infotypes a character vector containing the infotypes will be simulated in the additional bootstraps;
#'                  or a list with each element as a character vector, used for doing multiple additional
#'                  bootstraps with different settings of infotypes.
#' @param seed a numerical value indicates the seed for the randomness.
#'
#' @return
#' @export
#'
#' @examples
boot.acre = function(fit, N = 30, n.cores = 1, infotypes = NULL, seed = NULL){
  stopifnot(N > 2)
  fit.og = fit
  arguments = fit$args

  #deal with infotypes

  if(!is.null(infotypes)){
    #all possible infotypes
    all_infotypes = fit$infotypes

    if(is.null(all_infotypes)){
      warning('argument "infotypes" will be ignored as there is no extra infomation in the original model.')
      infotypes = NULL
    } else {
      #check it is a list or character vector
      stopifnot(any(is(infotypes, 'character'), is(infotypes, 'list')))
      if(is(infotypes, 'list')) stopifnot(all(sapply(infotypes, function(x) is(x, 'character'))))

      #if it is a character vector, convert it to a list to make future operation easier
      if(is(infotypes, 'character')) infotypes = list(infotypes)

      if(any(!do.call('c', infotypes) %in% all_infotypes)){
        stop(paste0('argument infotypes only accept following characters: ',
                    paste(all_infotypes, collapse = ', ')))
      }
    }
  }


  #since in parametric bootstrap we know the true values of each parameter
  #there are benefits to use the true values as start values
  sv.og = arguments$sv

  arguments$sv = get_sv_for_boot(fit)

  #extract some important components from original fit, maybe they need be adjusted
  #or they need be retained
  dims= get_dims_tmb(fit)
  traps = arguments$traps
  masks = arguments$mask



  ###############################################################################
  ##mrds, I'm not sure how to deal with 'mrds' in bootstrap yet, because the number
  ##of detection from simulation is random, however 'mrds' has provided the number
  ##of detection with their location. how to solve this conflict?
  #temporarily force it to be FALSE since not sure how to deal with it yet.
  if(fit$fit.types['mrds']){
    fit$infotypes = fit$infotypes[which(fit$infotypes!='mrds')]
    if(fit$n.sessions == 1){
      arguments$capt[['mrds']] = NULL
    } else {
      for(s in 1:fit$n.sessions) arguments$capt[[s]][['mrds']] = NULL
    }
  }

  if(!is.null(infotypes)){
    for(i in 1:length(infotypes)){
      if('mrds' %in% infotypes[[i]]) infotypes[[i]] = infotypes[[i]][which(infotypes[[i]] != 'mrds')]
    }
  }


  fit$fit.types['mrds'] = FALSE
  ##############################################################################
  is.mrds = fit$fit.types['mrds']
  if(is.mrds){
    if(fit$n.sessions == 1){
      mrds.loc = arguments$capt[['mrds']]
    } else {
      mrds.loc = lapply(arguments$capt, function(x) x[['mrds']])
    }
  }



  #coefficients names
  coefs = coef(fit, 'linked')
  n_pars = length(coefs)
  par_names = names(coefs)

  #set seed
  if(!is.null(seed)){
    set.seed(seed)
  } else {
    set.seed(sample(1:1e8, size = 1))
  }

  seed_boot = sample(1:1e8, size = N)


  if(n.cores == 1){
    ncol_res = n_pars + 1 + dims$n.sessions
    colnames_res = c(par_names, 'maxgrad', paste('esa', seq(dims$n.sessions), sep = "_"))


    ########################################################
    res = boot_step(seed = seed_boot,
                    N = N,
                    fit = fit,
                    arguments = arguments,
                    dims = dims,
                    infotypes = fit$infotypes,
                    len_output = ncol_res,
                    name_output = colnames_res)
    ########################################################


    #Additional bootstraps
    if(!is.null(infotypes)){
      extra.res <- vector(mode = "list", length = length(infotypes))
      names(extra.res) <- names(infotypes)

      for (i in seq(from = 1, by = 1, along.with = infotypes)){
        new_args <- arguments
        new_args$capt <- arguments$capt[c("bincapt", infotypes[[i]])]
        new_fit <- suppressWarnings(do.call("fit_og", new_args))
        tem = coef(new_fit, 'linked')
        new_ncol_res = length(tem) + 1 + dims$n.sessions
        new_colnames_res = c(names(tem), 'maxgrad')

        ####################################################################################
        extra.res[[i]] <- suppressWarnings(boot_step(seed = seed_boot,
                                                     N = N,
                                                     fit = new_fit,
                                                     arguments = new_args,
                                                     dims = dims,
                                                     infotypes = infotypes[[i]],
                                                     len_output = new_ncol_res,
                                                     name_output = new_colnames_res))
        ####################################################################################
        colnames(extra.res[[i]]) <- new_colnames_res

      }

      #end of if(!is.null(infotypes))
    } else {
      extra.res = NULL
    }

    #end of n.core == 1
  } else {
    stop('parallel computing is not avaiable yet.')
  }

  #"res" contains all covariates under "link" scale
  #browser()
  tem = res_split(res, dims$n.sessions)
  res= tem$res
  maxgrads = tem$maxgrads
  res_esa = tem$res_esa

  out <- fit.og
  boot <- list(boots = res, maxgrads = maxgrads, res_esa = res_esa, extra.boots = extra.res)
  out$boot <- boot
  class(out) <- c("acre_boot", class(fit))
  return(out)
}


boot_step = function(seed, N, fit, arguments, dims, infotypes, len_output, name_output){
  animal.model = "mu" %in% names(get_coef(fit))
  set.seed(seed)
  output = matrix(NA, nrow = N, ncol = len_output)
  #main bootstrap
  tem = suppressMessages(sim.capt(fit = fit, n.rand = N))
  capture_sim = tem$capt
  cue_rates_sim = tem$sim_cue_rates
  #browser()

  for(n in 1:N){
    if(nrow(capture_sim[[n]]) > 0){
      #browser()
      arguments$capt = get_capt_for_boot(capture_sim[[n]], dims, infotypes, animal.model)
      arguments$cue.rates = cue_rates_sim[[n]]
      fit_boot = suppressWarnings(try(do.call('fit_og', arguments), silent = TRUE))
      #If unconverged, refit model with default start values.
      if ("try-error" %in% class(fit_boot) || fit_boot$maxgrad < -0.01){
        arguments$sv <- NULL
        fit_boot <- suppressWarnings(try(do.call("fit_og", arguments), silent = TRUE))
      }
      #If still unconverged, give up and report NA, which is the default value of the output,
      #so we need to do nothing. And if it converged, assign result to the output
      if (!("try-error" %in% class(fit_boot) || fit_boot$maxgrad < -0.01)){
        output[n, ] = c(coef(fit_boot, 'linked'), fit_boot$maxgrad, coef(fit_boot, 'derived'))
      }

    } else {
      #default of output is NA, so no need to modify the columns of parameters excluding "D" related parameter
      #change the intercept or its equivalent to -Inf, and other "D" related coefficients to 0
      i_D_int = which(name_output %in% c('D_link', 'D.(Intercept)_link'))
      i_D_other = setdiff(which(grepl("^D", name_output)), i_D_int)
      output[n, i_D_int] = -Inf
      output[n, i_D_other] = 0
    }
  }

  colnames(output) = name_output

  return(output)
}





