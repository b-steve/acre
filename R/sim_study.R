#' Simulate acre-formatted SCR data.
#'
#' @param sim.name a string denoting type of data to be simulated. Supported 
#'                 types can be found using the `get.dataset.names()` 
#'                 function.
#' @param n.rand a numeric value denoting the number of data sets to be 
#'               generated. If n.rand > 1, `output$capt` will be a list of 
#'               length `n.rand`.
#' @param seed a numeric value denoting the random seed, defaults to 810.
#' @param suppress.messages a logical value which indicates whether to suppress 
#'                          helper messages.
#' @param proportion.missing a numeric value used to set a proportion of the 
#'                           covariate data generated to NA (currently only 
#'                           bearing, distance and toa supported)
#' @param fit an object generated from the model fitting function "fit.acre()" or
#'            the bootstrap process "boot.acre()". If `fit` is provided, then all
#'            parameters will be taken from the fitted object.
#'
#' @return A data frame containing the simulated data, as well as arguments used.
#' @export
sim.data <- function(sim.name, n.rand, fit=NULL, seed = 810, suppress.messages = FALSE, 
                    proportion.missing=0) {
  
  if (!sim.name %in% get.dataset.names() && is.null(fit)) {
    stop('invalid input for "sim.name", which should be one of the following: \n"', 
         paste(get.dataset.names(), collapse = '", "'), '".\n\n',
         "Alternatively provide a fitted model using the \"fit\" argument.")
  }
  if (!is.numeric(proportion.missing) || proportion.missing < 0 || proportion.missing > 1) {
    stop("proportion.missing must be a numeric value between 0 and 1.")
  }
  stopifnot(is.numeric(n.rand))
  
  set.seed(seed)
  
  if (!suppress.messages) 
    message("Simulation progress:")
  
  ## Generate the simulation arguments, and then simulate the raw capture data.
  set.seed(seed)
  
  ## If fit is provided, then use the fitted parameters; otherwise use the default ones.
  sim_args <- if(is.null(fit)) sim.args.generator(sim.name) else list(fit = fit)
  
  sim_args$n.rand <- n.rand
  simulated_capt <- do.call('sim.capt', sim_args)
  
  ## Generate acre-formatted data arguments.
  output <- fit.args.generator.from.sim(sim.name, simulated_capt$args)
  
  if(is.null(fit)) {
    output$param <- sim_args$param
  } else {
    ## Get the data frame that contains the link function for each parameter.
    dat_par <- get.data.param(fit)
    ## Then get the parameter values before back-transforming.
    param <- get.coef(fit)
    ## The "param" argument in the simulation requires back-transformed values.
    output$param <- param.transform(param, dat_par, back = TRUE)
  }
  
  if (!suppress.messages) 
    message("Converting raw data to acre format...")
  
  if (n.rand == 1) {
    ## Set data to NA if required.
    n_missing <- floor(nrow(simulated_capt$capt) * proportion.missing)
    if (n_missing > 0) {
      simulated_capt$capt <- set.detection.data.NA(simulated_capt$capt, 
                                                   sim.name, n_missing)
    }
    
    output$capt <- create.capt(simulated_capt$capt, output$traps)
  } else {
    ## Note that if we want multiple data sets, they will be stored as
    ## a list in output$capt.
    for (i in 1:n.rand) {
      ## Set data to NA if required.
      n_missing <- floor(nrow(simulated_capt$capt[[i]]) * proportion.missing)
      simulated_capt$capt[[i]] <- set.detection.data.NA(simulated_capt$capt[[i]], 
                                                   sim.name, n_missing)
      
      output$capt[[i]] <- create.capt(simulated_capt$capt[[i]], output$traps)
    }
  }
  
  if (!suppress.messages) 
    message(paste0("Successfully simulated ", n.rand, " capture data sets."))
  
  return(output)
}

#' Runs a simulation study 
#'
#' @param sim.name a string denoting type of data to be simulated. Supported 
#'                 types can be found using the `get.dataset.names()` 
#'                 function.
#' @param n.rand a numeric value denoting the number of data sets to be 
#'               generated and fit. Must be > 1.
#' @param n.cores a numeric value denoting the number of cores to use for 
#'                simulation. Defaults to `floor(parallel::detectCores() * 0.80)`
#' @param save a logical value, indicating whether the sim study data should be saved.
#'             By default creates a "sim_study" folder to save RDS object to. 
#' @param plot a logical value, indicating whether or not to plot the simulated 
#'             estimates.
#' @param proportion.missing a numeric value used to set a proportion of the
#'                           covariate data to NA. Will apply to all covariates 
#'                           present.
#' @param fit an object generated from the model fitting function "fit.acre()" or
#'            the bootstrap process "boot.acre()". If `fit` is provided, then all
#'            parameters will be taken from the fitted object.
#'
#' @return A list of length `n.rand` containing fitted coefficients from the 
#'         simulation study.
#' @export
sim.study <- function(sim.name, n.rand, fit=NULL, n.cores=NULL, save=T, plot=T, 
                     proportion.missing=0) {
  stopifnot(is.numeric(n.rand))
  stopifnot("simulation study must be run with at least 2 data sets. 
            If you wish to fit a single dataset consider using the `fit.acre()` 
            function" = n.rand > 1)
  if (!is.numeric(proportion.missing) || proportion.missing < 0 || proportion.missing > 1) {
    stop("proportion.missing must be a numeric value between 0 and 1.")
  }
  
  ## Simulate capture data.
  simulated_data <- sim.data(sim.name, n.rand, fit=fit, proportion.missing=proportion.missing)

  ## Generate the data frame that contains the default link function for each parameter.
  dat_par <- default.df.link()
  
  ## Get the linked values for all parameters.
  true_values <- param.transform(simulated_data$param, dat_par)
  
  ## Set up parallel computing.
  num_cores <- if(is.null(n.cores)) floor(parallel::detectCores() * 0.80) else n.cores
  message(paste0("Running simulation study using ", num_cores, " cores"))
  cl <- parallel::makeCluster(num_cores)
  
  est_values <- parallel::parLapply(cl, 1:n.rand, function(i) {
    library(acre)  ## Load the necessary package inside the cluster.
    
    fit_args <- simulated_data
    fit_args$capt <- fit_args$capt[[i]]
    fit_args$tracing <- FALSE
    
    sim_fit <- tryCatch({
      fit <- do.call('fit.from.settings', fit_args)
      fit$error <- FALSE
      fit
    }, error = function(e) {
      list(error = TRUE, message = e$message)
    })
    
    if (sim_fit$error) {
      return(list(error = TRUE, message = sim_fit$message))
    } else {
      est_value_linked <- cbind(Est = coef(sim_fit), 
                                Std = stdEr(sim_fit), confint(sim_fit))
      
      est_value_fitted <- cbind(Est = coef(sim_fit, types = "fitted"), 
                                Std = stdEr(sim_fit, types = "fitted"), 
                                confint(sim_fit, types = "fitted"))
      
      return(list(fitted = est_value_fitted, linked = est_value_linked))
    }
  })
  
  message("Finished sim study successfully.")
  parallel::stopCluster(cl)
  est_values[["true_link_parameters"]] <- true_values
  
  ## Create a "sim_study" folder if it does not exist.
  dir.create("sim_study", showWarnings = FALSE)
  
  if (save) {
    save_name <- paste0("sim_study/", sim.name, 
                        "_missing_", proportion.missing, "_n_", n.rand)
    
    saveRDS(est_values, save_name)
    message(paste("Saved study in", save_name))
  }
  
  if (plot) {
    ## Plot est_values.
    for(i in names(true_values)){
      tru <- true_values[[i]]
      i <- paste0(i, "_link")
      est <- sapply(est_values[1:(length(est_values)-1)], function(x) x[["linked"]][[i, 1]])
      

      hist(est, xlim = range(est, tru),
           main = paste0(i, ': Sim vs. True'),
           xlab = 'Value')

      abline(v = tru, col = 2)
      legend('topright', 'true value', col = 2, lty = 1)
      box()
      
    }
  }

  return(est_values)
}


## Set the first n additional detection data values to NA.
set.detection.data.NA <- function(capt.data, dataset.name, n.missing) {
  if (n.missing == 0) return (capt.data)
  
  detection_data_types <- c("bearing", "dist", "toa")
  
  ## The mul_ses data sets contain bearing and distance info, so just update the name
  ## to make sure they get appropriately set to NA.
  if (dataset.name %in% c("mul_ses", "mul_ses_ext")) {
    dataset.name <- paste0("bearing_dist_", dataset.name)
  }
  
  ## For each type of additional data.
  for (data_type in detection_data_types) {
    
    ## If the data contain this type of additional info.
    if (grepl(data_type, dataset.name)) {
      
      ## If all the data are NA, just remove the covariate column entirely.
      if (n.missing == nrow(capt.data)) {
        capt.data[data_type] <- NULL
      }
      else {
        ## Set the desired number of rows to NA.
        capt.data[sample(1:nrow(capt.data), n.missing), data_type] <- NA
      }
      
    }
  }
  
  return(capt.data)
}

