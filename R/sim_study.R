#' Simulate acre-formatted SCR data.
#'
#' @param sim_name a string denoting type of data to be simulated. Supported 
#'                 types can be found using the \code{get_dataset_names()} 
#'                 function.
#' @param n.rand a numeric value denoting the number of data sets to be 
#'               generated. If n.rand > 1, \code{output$capt} will be a list of 
#'               length \code{n.rand}.
#' @param seed a numeric value denoting the random seed, defaults to 810.
#' @param suppress_messages a logical value which indicates whether to suppress 
#'                          helper messages.
#' @param proportion_missing a numeric value used to set a proportion of the 
#'                           covariate data generated to NA (currently only 
#'                           bearing, distance and toa supported)
#' @param fit an object generated from the model fitting function "fit.acre_tmb()" or
#'            the bootstrap process "boot.acre()". If `fit` is provided, then all
#'            parameters will be taken from the fitted object.
#'
#' @return A data frame containing the simulated data, as well as arguments used.
#' @export
#'
#' @examples simulated_data <- sim_data("bearing_hn", 5)
sim_data = function(sim_name, n.rand, fit=NULL, seed = 810, suppress_messages = F, 
                    proportion_missing=0) {
  
  if (!sim_name %in% get_dataset_names() && is.null(fit)) {
    stop('invalid input for "sim_name", which should be one of the following: \n"', 
         paste(get_dataset_names(), collapse = '", "'), '".\n\n',
         "Alternatively provide a fitted model using the \"fit\" argument.")
  }
  if (!is.numeric(proportion_missing) || proportion_missing < 0 || proportion_missing > 1) {
    stop("proportion_missing must be a numeric value between 0 and 1.")
  }
  stopifnot(is.numeric(n.rand))
  
  set.seed(seed)
  
  if (!suppress_messages) 
    message("Simulation progress:")
  
  # Generate the simulation arguments, and then simulate the raw capture data
  set.seed(seed)
  
  # If fit is provided then use the fit parameters, otherwise use the default ones
  sim_args = if(is.null(fit)) sim_args_generator(sim_name) else list(fit = fit)
  
  sim_args$n.rand = n.rand
  simulated_capt = do.call('sim.capt', sim_args)
  
  # Generate acre formatted data arguments
  output = fit_args_generator_from_sim(sim_name, simulated_capt$args)
  
  if(is.null(fit)) {
    output$param = sim_args$param
  } else {
    # Get the data frame which contains link function for each parameter
    dat_par = get_data_param(fit)
    # Then get the parameter's values before back transforming
    param = get_coef(fit)
    # The "param" in simulation requires back transformed value
    output$param = param_transform(param, dat_par, back = TRUE)
  }
  
  if (!suppress_messages) 
    message("Converting raw data to acre format...")
  
  if (n.rand == 1) {
    # Set data to NA if required
    n_missing <- floor(nrow(simulated_capt$capt) * proportion_missing)
    if (n_missing > 0) {
      simulated_capt$capt <- set_detection_data_NA(simulated_capt$capt, 
                                                   sim_name, n_missing)
    }
    
    output$capt = create.capt(simulated_capt$capt, output$traps)
  } else {
    # Note that in the case we want multiple data sets, they will be stored as 
    # a list in output$capt
    for (i in 1:n.rand) {
      # Set data to NA if required
      n_missing <- floor(nrow(simulated_capt$capt[[i]]) * proportion_missing)
      simulated_capt$capt[[i]] <- set_detection_data_NA(simulated_capt$capt[[i]], 
                                                   sim_name, n_missing)
      
      output$capt[[i]] = create.capt(simulated_capt$capt[[i]], output$traps)
    }
  }
  
  if (!suppress_messages) 
    message(paste0("Successfully simulated ", n.rand, " capture data sets."))
  
  return(output)
}

#' Runs a simulation study 
#'
#' @param sim_name a string denoting type of data to be simulated. Supported 
#'                 types can be found using the \code{get_dataset_names()} 
#'                 function.
#' @param n.rand a numeric value denoting the number of data sets to be 
#'               generated and fit. Must be > 1.
#' @param n.cores a numeric value denoting the number of cores to use for 
#'                simulation. Defaults to \code{floor(parallel::detectCores() * 0.80)}
#' @param save a logical value, indicating whether the sim study data should be saved.
#'             By default creates a "sim_study" folder to save RDS object to. 
#' @param plot a logical value, indicating whether or not to plot the simulated 
#'             estimates.
#' @param proportion_missing a numeric value used to set a proportion of the
#'                           covariate data to NA. Will apply to all covariates 
#'                           present.
#' @param fit an object generated from the model fitting function "fit.acre_tmb()" or
#'            the bootstrap process "boot.acre()". If `fit` is provided, then all
#'            parameters will be taken from the fitted object.
#'
#' @return A list of length `n.rand` containing fitted coefficients from the 
#'         simulation study.
#' @export
#'
#' @examples study <- sim_study("bearing_hn", 5, save=F)
#' @examples study_from_fit <- sim_study("example_model", fit=example.fit, 5, save=F)
sim_study = function(sim_name, n.rand, fit=NULL, n.cores=NULL, save=T, plot=T, 
                     proportion_missing=0) {
  stopifnot(is.numeric(n.rand))
  stopifnot("simulation study must be run with at least 2 data sets. 
            If you wish to fit a single dataset consider using the `fit.acre()` 
            function" = n.rand > 1)
  if (!is.numeric(proportion_missing) || proportion_missing < 0 || proportion_missing > 1) {
    stop("proportion_missing must be a numeric value between 0 and 1.")
  }
  
  # Simulate capture data
  simulated_data <- sim_data(sim_name, n.rand, fit=fit, proportion_missing=proportion_missing)

  # generate the data frame which contains default link function for each parameter
  dat_par = default_df_link()
  
  # get the linked values for all parameters
  true_values = param_transform(simulated_data$param, dat_par)
  
  # Setup parallel computing stuff
  num_cores <- if(is.null(n.cores)) floor(parallel::detectCores() * 0.80) else n.cores
  message(paste0("Running simulation study using ", num_cores, " cores"))
  cl <- parallel::makeCluster(num_cores)
  
  est_values <- parallel::parLapply(cl, 1:n.rand, function(i) {
    library(acre)  # Load necessary package inside the cluster
    
    fit_args <- simulated_data
    fit_args$capt <- fit_args$capt[[i]]
    fit_args$tracing <- FALSE
    
    sim_fit <- tryCatch({
      fit <- do.call('fit_og', fit_args)
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
      
      est_value_fitted <- cbind(Est = coef(sim_fit, type = "fitted"), 
                                Std = stdEr(sim_fit, type = "fitted"), 
                                confint(sim_fit, type = "fitted"))
      
      return(list(fitted = est_value_fitted, linked = est_value_linked))
    }
  })
  
  message("Finished sim study successfully.")
  parallel::stopCluster(cl)
  est_values[["true_link_parameters"]] <- true_values
  
  # Create a "sim_study" folder if it does not exist
  dir.create("sim_study", showWarnings = FALSE)
  
  if (save) {
    save_name <- paste0("sim_study/", sim_name, 
                        "_missing_", proportion_missing, "_n_", n.rand)
    
    saveRDS(est_values, save_name)
    message(paste("Saved study in", save_name))
  }
  
  if (plot) {
    # plot the est_values
    for(i in names(true_values)){
      tru = true_values[[i]]
      i <- paste0(i, "_link")
      est = sapply(est_values[1:(length(est_values)-1)], function(x) x[["linked"]][[i, 1]])
      

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


#' Sets the first n additional detection data information to NA
#'
#' @param data Data set to be adjusted
#' @param dataset_name Name of the data set (must include the additional 
#'                    information to be removed in the name)
#' @param n_missing Number of data to be set to NA
#' @param is_sim True if it is simulation data
#'
#' @return 
#'
#' @examples
set_detection_data_NA <- function(capt_data, dataset_name, n_missing) {
  if (n_missing == 0) return (capt_data)
  
  detection_data_types <- c("bearing", "dist", "toa")
  
  # mul_ses data sets contain bearing and distance info, so just update name
  # to make sure they get appropriately set to NA
  if (dataset_name %in% c("mul_ses", "mul_ses_ext")) {
    dataset_name <- paste0("bearing_dist_", dataset_name)
  }
  
  # For each type of additional data
  for (data_type in detection_data_types) {
    
    # If data contains this type of additional info
    if (grepl(data_type, dataset_name)) {
      
      # If all the data is NA, just remove the covariate column entirely
      if (n_missing == nrow(capt_data)) {
        capt_data[data_type] <- NULL
      }
      else {
        # Set the desired number of rows to NA
        capt_data[sample(1:nrow(capt_data), n_missing), data_type] <- NA
      }
      
    }
  }
  
  return(capt_data)
}

#' Run a simulation study with missing covariate data.
#'
#' @param dataset_name string denoting the name of the dataset to be studied.
#' @param n.rand numeric value indicating the number of datasets to be simulated.
#'
#' @return 
#'
#' @examples sim_study_for_missing_data("bearing_hn", 20)
sim_study_for_missing_data <- function(dataset_name, n.rand) {
  # for (p_missing in seq(0,1,length=11)) {
  for (p_missing in seq(0,1,length=11)) {
    sim_study_updated(dataset_name, n.rand, T, proportion_missing = p_missing)
  }
}
