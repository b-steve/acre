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
#'
#' @return A data frame containing the simulated data, as well as arguments used.
#' @export
#'
#' @examples simulated_data <- sim_data("bearing_hn", 5)
sim_data = function(sim_name, n.rand, seed = 810, suppress_messages = F, 
                    proportion_missing=0){
  # Input validation
  if (!sim_name %in% get_dataset_names()) {
    stop('invalid input for "sim_name", which should be one of the following: "', 
         paste(get_dataset_names(), collapse = '", "'), '"')
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
  sim_args = sim_args_generator(sim_name)
  sim_args$n.rand = n.rand
  simulated_capt = do.call('sim.capt', sim_args)
  
  # Generate acre formatted data arguments
  output = fit_args_generator_from_sim(sim_name, simulated_capt$args)
  output$param = sim_args$param
  
  # Set data to NA if required

  
  if (!suppress_messages) 
    message("Converting raw data to acre format...")
  
  if (n.rand == 1) {
    n_missing <- floor(nrow(simulated_capt$capt) * proportion_missing)
    simulated_capt$capt <- set_detection_data_NA(simulated_capt$capt, 
                                                 sim_name, n_missing)
    
    output$capt = create.capt(simulated_capt$capt, output$traps)
  } else {
    # Note that in the case we want multiple data sets, they will be stored in 
    # a list in output$capt
    for (i in 1:n.rand) {
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

#' Updated simulation study
#'
#' @param sim_name 
#' @param n.rand 
#'
#' @return
#' @export
#'
#' @examples
sim_study_updated = function(sim_name, n.rand, save=T, proportion_missing=0) {
  stopifnot(is.numeric(n.rand))
  stopifnot("simulation study must be run with at least 2 data sets. 
            If you wish to fit a single dataset use the `fit.acre()` 
            function" = n.rand > 1)
  if (!is.numeric(proportion_missing) || proportion_missing < 0 || proportion_missing > 1) {
    stop("proportion_missing must be a numeric value between 0 and 1.")
  }
  
  # Simulate capture data
  simulated_data <- sim_data(sim_name, n.rand, proportion_missing=proportion_missing)

  # generate the data frame which contains default link function for each parameter
  dat_par = default_df_link()
  
  # get the linked values for all parameters
  true_values = param_transform(simulated_data$param, dat_par)
  
  # Setup parallel computing stuff
  num_cores <- floor(parallel::detectCores() * 0.5)
  message(paste0("Running simulation study using ", num_cores, " cores"))
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`
  
  est_values <- foreach::foreach(i = 1:n.rand, .packages = 'acre', 
                                 .errorhandling = 'pass') %dopar% {
 # for(i in 1:n.rand) {
                                   
    fit_args <- simulated_data
    fit_args$capt <- fit_args$capt[[i]]
    
    fit_args$tracing = FALSE

    sim_fit <- tryCatch({
      fit <- do.call('fit_og', fit_args)
      fit$error <- FALSE
      fit
    }, error = function(e) {
      # message(paste("Error in simulation", i, ":", e))
      list(error = TRUE, message = e$message)
    })
    
    # sim_fit <- do.call('fit_og', fit_args)
    
    if (sim_fit$error) {
      return(list(error = TRUE, message = sim_fit$message))
    } else {
      est_value_linked <- cbind(Est = coef(sim_fit), 
                                Std = stdEr(sim_fit), confint(sim_fit))
      
      est_value_fitted <- cbind(Est = coef(sim_fit, type="fitted"), 
                                Std = stdEr(sim_fit, type="fitted"), 
                                confint(sim_fit, type="fitted"))
      
      # Adds to the current list of est_values (doParallel trickery)
      return(list(fitted = est_value_fitted, linked = est_value_linked))
    }
  }
  
  message("Finished sim study successfully.")
  
  parallel::stopCluster(cl)
  
  est_values[["true_link_parameters"]] <- true_values
  
  if (save) {
    saveRDS(est_values, paste0("test_fits/", sim_name, 
                            "_missing_", proportion_missing, "_n_", n.rand))
  }

  return(est_values)
}

#' Simulates, fits and plots capture-recapture data
#'
#' @param sim_name Name of data set to be simulated
#' @param n.rand Number of data sets to simulate and fit
#' @param fit True if we wish to fit and plot the simulation data
#' @param seed Random seed
#' @param proportion_missing The proportion of additional detection data to be 
#'                            removed. Defaults to 0.
#'
#' @return
#' @export
#'
#' @examples
sim_study = function(sim_name, n.rand = 1, fit = FALSE, seed = 810, proportion_missing=0.0){
  sim_args = sim_args_generator(sim_name)
  sim_args$n.rand = n.rand
  set.seed(seed)
  if(n.rand != 1) message("simulation progress:")
  simulated_capt = do.call('sim.capt', sim_args)
  output = list()
  output$sim_capt = simulated_capt$capt
  output$sim_args = sim_args
  
  if(fit){
    message("[SIM_STUDY] Modelling started.")
    # sort out arguments required for model fitting first (excluding capt, as this could vary depends on n.rand)
    fit_args = fit_args_generator_from_sim(sim_name, simulated_capt$args)
    
    # if we only run 1 simulation and fit model to it, we could store the output of fit
    # if we run more simulations and fit models to all of them, we could not store that many
    # output objects, so we only plot the result of linked coefficients in a histogram with the
    # "true" value marked as a red vertical line
    if(n.rand == 1){
      fit_args$capt = create.capt(simulated_capt$capt, fit_args$traps)
      sim_fit = do.call('fit_og', fit_args)
      
      output$sim_fit_args = fit_args
      output$sim_fit = sim_fit
    } else {
      # generate the data frame which contains default link function for each parameter
      dat_par = default_df_link()
      
      # get the linked values for all parameters
      true_values = param_transform(sim_args$param, dat_par)
      est_values = vector('list', n.rand)
      fit_args$tracing = FALSE
      
      # Setup parallel computing stuff
      num_cores <- parallel::detectCores() - 1
      message(paste0("[SIM_STUDY] Running sim study using", num_cores, " cores"))
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      
      # for(i in 1:n.rand) {
      foreach::foreach(i = 1:n.rand, .packages = 'acre') %dopar% {
        # Removes data if necessary (none removed by default)
        n_missing <- floor(proportion_missing * nrow(simulated_capt$capt[[i]]))
        simulated_capt$capt[[i]] <- set_detection_data_NA(
          simulated_capt$capt[[i]],
          sim_name,
          n_missing,
          is_sim = T
        )
        
        # Fit model
        fit_args$capt = create.capt(simulated_capt$capt[[i]], fit_args$traps)
        sim_fit = do.call('fit_og', fit_args)
        est_values[[i]] = get_coef(sim_fit)
        message(paste0("[SIM_STUDY] finished: ", i, "/", n.rand))
        #write.csv(sim_fit$args$par.extend$data$mask, paste0('df_m_fit', i, '.csv'), row.names = F)
        
        if (save_fits) {
          save_sim_fit(sim_fit, i, sim_name, proportion_missing)
        }
      }
      
      parallel::stopCluster(cl)
      

      #remove the capture history and output the rest of arguments for model fitting
      fit_args$capt = NULL
      output$sim_fit_args = fit_args
      #we only output the linked scale of estimated coefficients
      output$sim_fit_coef_link = est_values


      #plot the est_values

      # for(i in names(true_values)){
      #   for(j in 1:length(true_values[[i]])){
      #     #extract estimations as a vector and the true value
      #     est = sapply(est_values, function(x) x[[i]][j])
      #     tru = true_values[[i]][j]
      # 
      #     hist(est, xlim = range(est, tru),
      #          main = paste0(i, '[', j, ']_link: Sim vs. True'),
      #          xlab = 'Value')
      # 
      #     abline(v = tru, col = 2)
      #     legend('topright', 'true value', col = 2, lty = 1)
      #     box()
      #   }
      # }
      
    }
  }
  
  return(output)
}

#' Saves models fitted to simulation data
#'
#' @param sim_fit The fitted model
#' @param i The iteration / simulation number
#' @param sim_name The simulation / data set name.
#' @param proportion_missing The proportion of additional detection data to be 
#'                            removed.
#'
save_sim_fit <- function(sim_fit, i, sim_name, proportion_missing=0) {
  # Save the data to test_fits/sim/sim_name_i
  save_to_location <- paste0("test_fits/sim/", sim_name, "/missing_",
                             proportion_missing, "/")
  save_to_name <- paste0(sim_name, "_", i)
  save_to <- paste0(save_to_location, save_to_name)

  conf_int <- confint.acre_tmb(sim_fit)
  coefs <- coef.acre_tmb(sim_fit)
  errs <- stdEr.acre_tmb(sim_fit)

  fit_data <- data.frame(conf_int = conf_int,
                         coefs = as.vector(coefs),
                         errs = as.vector(errs))

  # Update column names
  names(fit_data)[1] <- "conf_int_2.5"
  names(fit_data)[2] <- "conf_int_97.5"

  saveRDS(fit_data, save_to)
  
  message(paste0("Saved sim fit data to: ", save_to_location))
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
      # Set the desired number of rows to NA
      capt_data[sample(1:nrow(capt_data), n_missing), data_type] <- NA
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
#' @export
#'
#' @examples sim_study_for_missing_data("bearing_hn", 20)
sim_study_for_missing_data <- function(dataset_name, n.rand) {
  for (p_missing in seq(0,1,length=11)) {
    sim_study_updated(dataset_name, n.rand, T, proportion_missing = p_missing)
  }
}
