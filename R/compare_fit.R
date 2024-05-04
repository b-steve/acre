fit_models_with_missing_data <- function(dataset_names=NULL, remove_every_n) {
  # If no data-set names provided, fit all the test data-sets 
  # (also removing .rda extension)
  if (is.null(dataset_names)) {
    dataset_names <- gsub("\\.rda$", "", list.files("data"))
  }
  
  for (dataset_name in dataset_names) {
    data <- get(dataset_name)
    
    # Remove more data every iteration
    for (n_missing in seq(0, nrow(data$capt), remove_every_n)) {
      adj_data <- set_detection_data_NA(data, dataset_name, n_missing)
      
      # Convert data to acre readable format
      acre_adj_data <- read.acre(adj_data$capt, adj_data$traps, 
                control_create_mask = adj_data$control_create_mask)
      
      # Notice data$ss.opts will be NULL if not a ss model
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

compare_missing_data_via_sim_study <- function(n.rand, datasets) {
  for (dataset in datasets) {
    for (proportion_missing in seq(from=0.0, to=1, by=0.1)) {
      sim_study(dataset, 
                fit=T, 
                n.rand = n.rand, 
                proportion_missing = proportion_missing)
    }
  }
}

plot_bootstrap_fits <- function(folders=NULL) {
  base_folder <- "test_fits/sim/dist_hn/"  
  
  folders <- c("missing_0",
               "missing_0.1",
               "missing_0.2",
               "missing_0.3",
               "missing_0.4",
               "missing_0.5")
  
  
  coefs <- list()
  for (folder in folders) {
    path_name <- paste0(base_folder, folder)
    coefs_vector <- numeric()
    
    print(length(list.files(path_name)))
    
    for (file in list.files(path_name)[1:97]) {
      data <- readRDS(paste0(path_name, "/", file))
  
      coefs_vector <- c(coefs_vector, exp(data[2, "coefs"]))
    }
    
    coefs[[folder]] <- coefs_vector
  }
  
  # return (coefs)
  # plot(coefs_vector, rep(1, length(coefs_vector)), pch = 19, xlab = "Coefficient Value", ylab = "", main = "Plot of Coefficients")
  # Create a data frame for plotting
  plot_data <- data.frame(
    Coefficient = unlist(coefs),
    missing = (0:(length(folders)-1)) * 0.1
  )
  
  # Plotting
  plot(missing ~ Coefficient, data = plot_data, pch = 19,
       xlab = "Coefficient value", ylab = "Proportion missing", 
       main = "Plot of Coefficients")
  abline(v=5.33)
  
}






























