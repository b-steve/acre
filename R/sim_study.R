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

      for(i in names(true_values)){
        for(j in 1:length(true_values[[i]])){
          #extract estimations as a vector and the true value
          est = sapply(est_values, function(x) x[[i]][j])
          tru = true_values[[i]][j]

          hist(est, xlim = range(est, tru),
               main = paste0(i, '[', j, ']_link: Sim vs. True'),
               xlab = 'Value')

          abline(v = tru, col = 2)
          legend('topright', 'true value', col = 2, lty = 1)
          box()
        }
      }
      
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
set_detection_data_NA <- function(data, dataset_name, n_missing, is_sim = F) {
  if (n_missing == 0) return (data)
  
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
      
      # Note for sim data, we are passing the capture info directly
      if (is_sim) {
        # Set the desired number of rows to NA
        capt_data[1:n_missing, data_type] <- NA
      } else {
        data$capt[1:n_missing, data_type] <- NA
      }
    }
  }
  
  return(data)
}
