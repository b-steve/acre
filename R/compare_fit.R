fit_models_with_missing_data <- function(dataset_names=NULL) {
  # Fitting for each test data-set (also removing .rda extension)
  if (is.null(dataset_names)) {
    dataset_names <- gsub("\\.rda$", "", list.files("data"))
  }
  
  for (dataset_name in dataset_names) {
    data <- get(dataset_name)
    
    # Remove an additional row each loop
    for (n_missing in seq(0, nrow(data$capt), 5)) {
      adj_data <- set_detection_data_NA(data, dataset_name, n_missing)
      
      # Convert data to acre readable format
      acre_adj_data <- read.acre(adj_data$capt, adj_data$traps, 
                control_create_mask = adj_data$control_create_mask)
      
      # Notice data$ss.opts will be NULL if not ss model
      fit <- fit.acre(acre_adj_data, ss.opts = data$ss.opts)
      
      # Extract fit info
      conf_int <- confint.acre_tmb(fit)
      coefs <- coef.acre_tmb(fit)
      errs <- stdEr.acre_tmb(fit)
      
      # Save it to the respective folder
      fit_data <- data.frame(conf_int = conf_int, 
                             coefs = as.vector(coefs), 
                             errs = as.vector(errs))
      
      # Update column names
      names(fit_data)[1] <- "conf_int_2.5"
      names(fit_data)[2] <- "conf_int_97.5"
      
      fit_name <- paste0("test_fits/missing/", dataset_name, "/missing_", n_missing)
      print(fit_name)
      saveRDS(fit_data, fit_name)
    }
  }
}


set_detection_data_NA <- function(data, dataset_name, n_missing) {
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
      
      # Set the desired number of rows to NA
      data$capt[1:n_missing, data_type] <- NA
    }
  }
  
  return(data)
}

simulate_comparison_datasets <- function() {
  # Record true parameter values
  
  # For each level of data removal
  sim.capt
  # Simulate 100 data sets
  
  # Fit 100 models
  
  
}

remove_capt_data_from_sim <- function(capt_data, sim_name, n_missing) {
  if (n_missing == 0) return (capt_data)
  
  detection_data_types <- c("bearing", "dist", "toa")
  
  # mul_ses data sets contain bearing and distance info, so just update name
  # to make sure they get appropriately set to NA
  if (sim_name %in% c("mul_ses", "mul_ses_ext")) {
    sim_name <- paste0("bearing_dist_", dataset_name)
  }
  
  # For each type of additional data
  for (data_type in detection_data_types) {
    
    # If data contains this type of additional info
    if (grepl(data_type, sim_name)) {
      
      # Set the desired number of rows to NA
      capt_data[1:n_missing, data_type] <- NA
    }
  }
  
  return(capt_data)
}


compare_via_sim_study <- function(n.rand, datasets) {
  
  for (dataset in datasets) {
    for (proportion_missing in seq(from=0.0, to=1, by=0.1)) {
      sim_study(dataset, 
                fit=T, 
                n.rand = n.rand, 
                proportion_missing = proportion_missing)
    }
  }
}

plot_bootstrap_fits <- function() {
  base_folder <- "test_fits/sim/dist_hn/"  
  folders <- c("missing_0",
               "missing_0.1",
               "missing_0.2",
               "missing_0.3",
               "missing_0.4")
  
  
  coefs <- list()
  for (folder in folders) {
    path_name <- paste0(base_folder, folder)
    coefs_vector <- numeric()
    
    for (file in list.files(path_name)) {
      data <- readRDS(paste0(path_name, "/", file))
  
      coefs_vector <- c(coefs_vector, data[2, "coefs"])
    }
    
    coefs[[folder]] <- coefs_vector
  }
  
  # return (coefs)
  # plot(coefs_vector, rep(1, length(coefs_vector)), pch = 19, xlab = "Coefficient Value", ylab = "", main = "Plot of Coefficients")
  # Create a data frame for plotting
  plot_data <- data.frame(
    Coefficient = unlist(coefs),
    Folder = 1:length(folders)
  )
  
  # Plotting
  plot(Folder ~ Coefficient, data = plot_data, pch = 19,
       xlab = "Folder", ylab = "Coefficient Value", 
       main = "Plot of Coefficients")
  
  }





