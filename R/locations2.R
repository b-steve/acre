# plot_location_density(ind_bd.fit, animal_id=1, call_id=1, plot_types = c("combined"), 
# mask=secr::make.mask(ind_bearing_dist$traps, 15, poly=data.frame(x=c(-15, 20, 20, -15), y=c(25, 25, -15, -15))))

# plot_location_density(ind_bd.fit, animal_id=1, call_id=1, plot_types = c("combined"), 
# mask=secr::make.mask(ind_bearing_dist$traps, 15, poly=data.frame(x=c(-30, 40, 40, -30), y=c(35, 35, -30, -30))))

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#' @param density_types 
#' @param combine 
#' @param warn_toa 
#' @param new_cov_data 
#'
#' @return
#' @export
#'
#' @examples
location_density <- function(fit, call_id, animal_id=NULL, session=1, 
                              mask=get_mask(fit)[[session]], 
                              density_types=c(fit$infotypes,"capt"), 
                              combine=T, warn_toa=F, new_cov_data=NULL) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  # Note each density function returns a df with columns:
  #   ID, x, y, [density_name]_density
  density_list <- lapply(density_types, function(d_type) {
    switch(
      d_type,
      "pois" = pois_density(fit, call_id, session, mask, new_cov_data),
      "capt" = capt_density(fit, call_id, animal_id, session, mask),
      "toa"  = toa_density(fit, call_id, animal_id, session, mask,
                           warn=warn_toa),
      "dist" = dist_density(fit, call_id, animal_id, session, mask),
      "ss"   = ss_density(fit, call_id, animal_id, session, mask),
      "bearing" = bearing_density(fit, call_id, animal_id, session, mask),
      stop(paste0("Unknown density type: '", d_type, 
                  "', currently supported types are: 'capt', 'ss', 'dist',
                   'bearing', 'pois'."))
    )
  })
  
  # For clarity
  names(density_list) <- density_types
  
  # Merge all the data frames together
  if (length(density_list) == 0) {
    # Edge case - no density dataframes
    warning("No densities were calculated. Returning an empty dataframe.")
    return(data.frame())
  } else {
    merged_densitys <- density_list[[1]]
    
    # Merge each subsequent data frame
    for (i in seq_along(density_list)[-1]) {
      merged_densitys <- merge(
        x  = merged_densitys,
        y  = density_list[[i]],
        by = c("ID", "x", "y"), # Make sure to match up appropriate x, y & ID's
        all = TRUE  # Full outer join (keeps any rows which don't "match" up), 
                    # shouldn't matter though as the df's have same structure  
      )
    }
    
    if (combine) {
      if (ncol(merged_densitys) == 4) {
        warning("Only 1 density type calculated. 
                Combined density will be identical to this density column.")
      }
      # Create the combined density column by multiplying all density columns 
      merged_densitys$combined_density <- Reduce(`*`, merged_densitys[-c(1:3)])
    }
  
    # Return the final data frame
    return(merged_densitys)
  }
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#' @param new_cov_data 
#' @param plot_types 
#' @param joint 
#' @param nlevels 
#' @param ann 
#'
#' @return
#' @export
#'
#' @examples
plot_location_density <- function(fit, call_id="all", animal_id=NULL, session=1, 
                      mask=get_mask(fit)[[session]],  new_cov_data=NULL,
                      plot_types="combined", joint=F, nlevels=8, ann=T) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  # If we are using the default call_id "all"
  if (!is.numeric(call_id)) {
    call_id <- unique(get_capt_by_id(fit, call_id, animal_id, session)$ID)
  }
  
  # isTRUE() for clarity? Never used it before, but might start using it.
  # Check if we are using a mask different than the one used to fit the model.
  # If so, we need to make sure the new loc covariates are provided if IHD.
  is.new.mask <- !missing(mask) && !isTRUE(all.equal(mask, get_mask(fit)[[session]]))
  if (fit$fit.ihd & is.new.mask & is.null(new_cov_data)){
    stop("Covariate values for the new mask object must be provided via the 
         'newdata' argument.")
  }
  
  # Figure out which densities we want
  if (length(plot_types) == 1 && plot_types == "combined") {
    density_types <- c(fit$infotypes, "capt", "pois")
  } else {
    density_types <- plot_types[plot_types != "combined"]
  }
  
  # if (length(density_types) == 1 & 
  #     "combined" %in% plot_types & 
  #     length(fit$infotypes) > 1) {
  #   stop("'plot_types' = 'combined' is only possible when 
  #        at least two other density types are to be plotted.")
  # }
  
  # For each call
  for (i in seq_along(call_id)) {
    cid <- call_id[i]
    
    # Calculate combined density by default, rather than bloating function args
    call_densities <- location_density(fit, cid, animal_id, session, mask = mask,
                                       density_types = density_types,
                                       combine = "combined" %in% plot_types,
                                       warn_toa = F, 
                                       new_cov_data = new_cov_data)
    
    # 2) Initialize a base ggplot
    p <- ggplot(call_densities, aes(x = x, y = y)) +
      coord_sf(xlim = range(mask[,1]), ylim = range(mask[,2]))
    
    if (ann) {
      p <- p + ggtitle(paste("Location density -", substitute(fit))) +
        labs(subtitle = create_locations_plot_subtitle(cid, animal_id, session))
    }
    
    for (d_type in plot_types) {
      density_column <- paste0(d_type, "_density")
      
      # Skip if that column doesn't exist 
      #   (in case location_density() not working properly)
      if (!density_column %in% names(call_densities)) {
        warning("No column ", density_column, " in call_densities. Skipping.")
        next
      }
      
      contour_color <- switch(
        d_type,
        "capt" = "black",
        "bearing" = "red",
        "dist" = "green",
        "combined" = "blue",
        "ss" = "orange",
        "toa" = "purple",
        "gray50" # default
      )
      
      # Note that we can't calculate / plot TOA densities in the case of only a
      # single trap detection. If this is the case, toa_density() returns an 
      # equal density across mask, so just check for that, rather than reading
      # in all the capt data.
      if (d_type == "toa" & 
          length(unique(call_densities[[density_column]])) == 1) {
        # If we are ONLY plotting toa for 1 call stop, otherwise just warn 
        if (length(call_id) == 1 & length(plot_types) == 1) {
          stop(paste0("Could not plot toa density for call: ", cid, 
                      ". Calls must be detected on at least 2 traps to use toa data."))
        } else {
          warning("Could not plot toa density for call: ", cid, 
                  ". Calls must be detected on at least 2 traps to use toa data.")
          next
        }
      }
      
      # Add contours
      # If we are only plotting a single density type, use a filled contour, 
      # otherwise use line contours.
      p <- p +
        geom_contour_helper(call_densities[[density_column]], 
                            color = contour_color, 
                            filled = length(plot_types) == 1, 
                            nlevels = nlevels)
    }
    
    # Add traps
    p <- p + plot_traps(fit, cid, animal_id, session, circle_traps = T)
    
    # Add bearing arrows
    if ("bearing" %in% density_types) {
      p <- p + plot_bearing_arrows(fit, cid, animal_id, session)
    }
    
    # Add distance circles
    if ("dist" %in% density_types) {
      p <- p + plot_distance_circles(fit, cid, animal_id, session)
    }
    
    # Add toa indicators
    if ("toa" %in% density_types) {
      p <- p + plot_toa_order(fit, cid, animal_id, session)
    }
    
    # Add ss indicators (not sure how I feel about this one right now)
    # if ("ss" %in% density_types) {
    #   p <- p + plot_ss_indicators(fit, cid, animal_id, session)
    # }

    p <- p + ggplot2::theme_bw()
    
    # Remove legend and axis if necessary
    if (!ann) {
      p <- p + ggplot2::theme(legend.position="none",      
                              axis.text.x=ggplot2::element_blank(),
                              axis.text.y=ggplot2::element_blank(),
                              axis.title.x=ggplot2::element_blank(),
                              axis.title.y=ggplot2::element_blank())
    }
    
    # Display plot
    print(p)
    
    # Make sure to only prompt user if we are plotting more than one plot, 
    # and it is not the last plot in the list
    if (length(call_id) > 1 && cid != call_id[[length(call_id)]]) {
      prompt_user_for_next_plot()
    }
  }
  
  invisible(p)
}

#' Title
#'
#' @param dens 
#' @param probability_breaks 
#'
#' @return
#' @export
#'
#' @examples
get_density_breaks_from_probability <- function(dens, probability_breaks) {
  dens.sort <- sort(dens, decreasing = TRUE)
  probs.sort <- cumsum(dens.sort)/sum(dens.sort)
  prob.levels <- probability_breaks
  levels <- numeric(length(prob.levels))
  for (i in 1:length(prob.levels)){
    levels[i] <- dens.sort[which(abs(probs.sort - prob.levels[i]) ==
                                min(abs(probs.sort - prob.levels[i])))[1]]
  }
}

#' Title
#'
#' @param dens 
#' @param color 
#' @param filled 
#'
#' @return
#' @export
#'
#' @examples
geom_contour_helper <- function(dens, color, filled=T, nlevels=8, prob=F) {
  # Tried using pretty() but it inconsistent output is annoying
  dens_range <- range(dens)
  levels <- seq(from = dens_range[1], to = dens_range[2], length.out = nlevels+1)
  
  # It might seem strange constructing the scale arguments function like this,
  # but doing so allows us to leave the "labels" field blank in the case of 
  # non-probability contours (prob=F), and use the ggplot default legend labels.
  scale_args <- list()
  
  # If we want probability labels
  if (prob) {
    probability_labels <- character(nlevels)
    for (i in 1:nlevels){
      probability_labels[i] <- format(
        round(sum(dens[dens > levels[i]]) / sum(dens), 2), 
        nsmall = 2)
    }
    
    scale_args$name <- "Probability"
    scale_args$labels <- probability_labels
  } else {
    # Otherwise use the default density labels
    scale_args$name <- "Density"
  }
    
  if (filled) {
    # base_colors <- rep(color, nlevels)
    base_colors <- viridisLite::viridis(nlevels)
    alpha_values <- seq(0, 1, length.out = nlevels)
    faded_colors <- unlist(
      Map(function(col, a) scales::alpha(col, a), base_colors, alpha_values)
    )
    scale_args$values = unname(faded_colors)
    
    return(
      list(
        ggplot2::geom_contour_filled(
          aes(z = dens),
          colour = "white",
          # colour = NA,
          # bins = 8,
          breaks = levels,
          linewidth=0.5
        ),
        do.call(ggplot2::scale_fill_manual, scale_args)
      )
    )
  } else {
    return(ggplot2::geom_contour(
      aes(z = dens),
      colour = color,
      alpha = 0.7,
      breaks = levels,
      bins = 15
    ))
  }
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param circle_traps 
#'
#' @return
#' @export
#'
#' @examples
plot_traps <- function(fit, call_id, animal_id, session=1, circle_traps=T) {
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting traps is only supported for single calls."))
  }
  
  traps <- get_trap(fit)[[session]]
  
  trap_points <- geom_point(data = traps, mapping = aes(x = x, y = y), 
             size = 3, col="red", stroke=1, shape = 4)
  
  if (circle_traps) {
    # Grab appropriate capture data
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt = T)[,-1]
    
    trap_circles <- geom_point(
      data = subset(traps, bincapt == 1), aes(x = x, y = y),
      shape = 21, color = "red", size = 4.5, stroke = 1)
    
    return(list(trap_points, trap_circles))
  } else {
    return(trap_points)
  }
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#'
#' @return
#' @export
#'
#' @examples
plot_toa_order <- function(fit, call_id, animal_id=NULL, session=1) {
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting toa indicators is only supported for single calls."))
  }
  
  if (!("toa" %in% fit$infotypes)) {
    stop("'fit' does not contain toa information.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  
  # Need at least 2 trap detections to use toa
  if (nrow(capt) < 2) {
    return()
  } else{
    traps <- get_trap(fit)[[session]]
    toa_order <- rank(capt$toa, ties.method = "min")
    
    return(ggplot2::geom_text(data = traps[capt$trap,], 
                              aes(x = x, y = y, label = toa_order),
                              size = 2, vjust = 1.5, hjust=-2, colour = "black"))
  }
  

}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#'
#' @return
#' @export
#'
#' @examples
plot_bearing_arrows <- function(fit, call_id, animal_id=NULL, session=1) {
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting ebaring arrows is only supported for single calls."))
  }
  
  # Make sure call_id is a single numeric
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting bearing arrows is only supported for single calls."))
  }
  
  if (!("bearing" %in% fit$infotypes)) {
    stop("'fit' does not contain bearing information.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)[,-1]
  
  is.dist <- "dist" %in% fit$infotypes
  traps <- get_trap(fit)[[session]]
  activated_traps <- subset(traps, bincapt == 1)
  bearing <- capt$bearing
  dist <- capt$dist
  
  # Set appropriate arrow length. Default to 1/3 buffer radius.
  arrow_len <- rep(0.333 * get_buffer(fit)[session], length(capt$bearing))
  if (is.dist) {
    # Account for the possibility that dist may be NA for some captures
    arrow_len[!is.na(capt$dist)] <- capt$dist[!is.na(capt$dist)]
  }
  
  sinb <- sin(bearing)*arrow_len
  cosb <- cos(bearing)*arrow_len
  
  arrows.df <- data.frame(
    x = activated_traps[, 1],
    y = activated_traps[, 2],
    xend = activated_traps[, 1] + sinb,
    yend = activated_traps[, 2] + cosb
  )
  
  # Remove any NA bearing captures
  arrows.df <- arrows.df[complete.cases(arrows.df), ]
  
  return(geom_segment(data = arrows.df, 
                      mapping = aes(x = x, y = y, xend = xend, yend = yend),
                      col="red",
                      arrow = arrow(length = unit(0.2, "cm")) # Arrow head 
                      )
         ) 
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#'
#' @return
#' @export
#'
#' @examples
plot_distance_circles <- function(fit, call_id, animal_id=NULL, session=1) {
  # Make sure call_id is a single numeric
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting distance circles is only supported for single calls."))
  }
  
  if (!("dist" %in% fit$infotypes)) {
    stop("'fit' does not contain distance information.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)[,-1]
  dist <- capt$dist
  traps <- get_trap(fit)[[session]]
  activated_traps <- subset(traps, bincapt == 1)
  
  circles.df <- data.frame(
    x = activated_traps[, 1], y = activated_traps[, 2], size=dist
  )
  
  # Remove any NA distance captures
  circles.df <- circles.df[complete.cases(circles.df), ]
  
  return(geom_point(
    data = circles.df,
    aes(x = x, y = y),
    shape = 21,
    color = "red",
    size = circles.df$size,  
    stroke = 1 
  ))
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#'
#' @return
#' @export
#'
#' @examples
plot_ss_indicators <- function(fit, call_id, animal_id, session=1) {
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting ss indicators is only supported for single calls."))
  }
  
  if (!("ss" %in% fit$infotypes)) {
    stop("'ss' does not contain toa information, cannot plot ss indicators.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  traps <- get_trap(fit)[[session]]
  
  return(list(geom_point(data = traps[capt$trap,], 
                         mapping = aes(x = x, y = y, colour = capt$ss), 
                         size = 5), 
              guides(colour = guide_colourbar(order = 1, 
                                              title = "Signal Strength"))
  )
  )
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#' @param dists 
#'
#' @return
#' @export
#'
#' @examples
dist_density <- function(fit, call_id, animal_id=NULL, session=1, 
             mask=get_mask(fit)[[session]], dists=NULL) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  if (!("dist" %in% fit$infotypes)) {
    stop(paste0("'fit' does not contain 'dist' information. ",
                "Infotypes available for the fit provided are: ", 
                paste0(fit$infotypes, collapse=", ")))
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
  
  # Perhaps need to add some validation for dists & mask
  if (is.null(dists)) {
    traps <- get_trap(fit)[[session]]
    dists <- distances(traps, as.matrix(mask))
  }

  # Grab our estimated alpha parameter
  alpha <- coef(fit, type="fitted")[["alpha"]]
  betas <- alpha/dists
  
  # A list to store each call density, will be combined into one data 
  density_list <- vector("list", length(call_id))
  
  for (i in seq_along(call_id)) {
    cid <- call_id[i]
    
    # Grab the capt & bincapt for this call ID
    sub_capt <- subset(capt, ID == cid)
    # Make sure to remove ID column when doing calculations
    sub_bin  <- subset(bincapt, ID == cid)[,-1] 
    dist.capt <- sub_capt$dist
    
    # Initialize density matrix for this call
    mask.dens <- matrix(1, nrow = sum(sub_bin), ncol = nrow(mask))
    
    # Calculate density for each triggered trap
    for (j in 1:sum(sub_bin)) {
      if (!is.na(dist.capt[j])) {
        mask.dens[j, ] <- dgamma(dist.capt[j], shape = alpha, rate = betas[j, ])
      }
    }
    
    mask.dens <- apply(mask.dens, 2, prod)

    density_list[[i]] <- data.frame(
      ID = cid, x = mask[, 1], y = mask[, 2], dist_density = mask.dens
    )
  }

  return(do.call(rbind, density_list))
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#'
#' @return
#' @export
#'
#' @examples
bearing_density <- function(fit, call_id, animal_id=NULL, session=1, 
                         mask=get_mask(fit)[[session]]) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  if (!("bearing" %in% fit$infotypes)) {
    stop(paste0("'fit' does not contain 'bearing' information. ",
                "Infotypes available for the fit provided are: ", 
                paste0(fit$infotypes, collapse=", ")))
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
  traps <- get_trap(fit)[[session]]
  
  # Grab our estimated kappa parameter
  kappa <- coef(fit, type="fitted")[["kappa"]]
  
  # A list to store each call density, will be combined into one data 
  density_list <- vector("list", length(call_id))
  
  for (i in seq_along(call_id)) {
    cid <- call_id[i]
    
    # Grab the capt & bincapt for this call ID
    sub_capt <- subset(capt, ID == cid)
    # Make sure to remove ID column when doing calculations
    sub_bin  <- subset(bincapt, ID == cid)[,-1] 
    bearing.capt <- sub_capt$bearing
    
    # Initialize density matrix for this call
    mask.dens <- matrix(1, nrow = sum(sub_bin), ncol = nrow(mask))
    
    # Calculate bearings for each triggered trap
    mask.bearings <- bearings(subset(traps, sub_bin == 1), as.matrix(mask))
    
    # Calculate the density for each triggered trap
    for (j in 1:sum(sub_bin)) {
      # Make sure bearing isn't missing
      if (!is.na(bearing.capt[j])) {
        mask.dens[j, ] <- CircStats::dvm(bearing.capt[j], 
                                         mu = mask.bearings[j, ], kappa = kappa)
      }
    }
    
    mask.dens <- apply(mask.dens, 2, prod)
    
    density_list[[i]] <- data.frame(
      ID = cid, x = mask[, 1], y = mask[, 2], bearing_density = mask.dens
    )
  }
  
  return(do.call(rbind, density_list))
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#' @param dists 
#'
#' @return
#' @export
#'
#' @examples
toa_density <- function(fit, call_id, animal_id=NULL, session=1, 
                        mask=get_mask(fit)[[session]], dists=NULL, warn=T) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  if (!("toa" %in% fit$infotypes)) {
    stop(paste0("'fit' does not contain 'toa' information. ",
                "Infotypes available for the fit provided are: ", 
                paste0(fit$infotypes, collapse=", ")))
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
  traps <- get_trap(fit)[[session]]
  
  # Perhaps need to add some validation for dists & mask
  if (is.null(dists)) {
    traps <- get_trap(fit)[[session]]
    dists <- distances(traps, as.matrix(mask))
  }
  
  # Grab our estimated sigma.toa parameter
  sigma.toa <- coef(fit, type="fitted")[["sigma.toa"]]
  
  # A list to store each call density, will be combined into one data 
  density_list <- vector("list", length(call_id))
  
  # Set to true if there was a capture history with only 1 detection
  # (need at least 2 detections to be able to use toa)
  single_trap_detection = F
  for (i in seq_along(call_id)) {
    cid <- call_id[i]
    
    # Grab the capt & bincapt for this call ID
    sub_capt <- subset(capt, ID == cid)
    # Make sure to remove ID column when doing calculations
    sub_bin  <- subset(bincapt, ID == cid)[,-1] 
    toa.capt <- sub_capt$toa
    
    if (sum(sub_bin) < 2) {
      single_trap_detection = T
      mask.dens <- 1
    } else {
      dists.mask <- dists[sub_bin == 1, ]
      prod.times <- toa.capt - dists.mask/fit$args$sound.speed
      toa.ssq <- apply(prod.times, 2, function(x) sum((x - mean(x))^2))
      mask.dens <- (2*pi*sigma.toa^2)^((1 - sum(sub_bin))/2)*
        exp(toa.ssq/(-2*sigma.toa^2))
    }

    density_list[[i]] <- data.frame(
      ID = cid, x = mask[, 1], y = mask[, 2], toa_density = mask.dens
    )
  }
  
  if (warn & single_trap_detection) {
    warning(paste("Some calls were detected by only a single trap. 
                  The use of toa data requires a call to be detected by at least 
                  2 traps. Toa density for corresponding calls has beem set to 1 
                  across the entire mask for such cases."))
  }
  return(do.call(rbind, density_list))
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#' @param dists 
#'
#' @return
#' @export
#'
#' @examples
ss_density <- function(fit, call_id, animal_id=NULL, session=1, 
                        mask=get_mask(fit)[[session]], dists=NULL) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  if (!("ss" %in% fit$infotypes)) {
    stop(paste0("'fit' does not contain 'ss' information. ",
                "Infotypes available for the fit provided are: ", 
                paste0(fit$infotypes, collapse=", ")))
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
  
  # Grab our estimated parameters, along with other model fitting parameters
  cutoff <- get_ss.opts(fit)$cutoff
  sigma.ss <- coef(fit, type="fitted")[["sigma.ss"]]
  b1.ss <- coef(fit, type="fitted")[["b1.ss"]]
  b0.ss <- coef(fit, type="fitted")[["b0.ss"]]
  detfn <- get_detfn(fit)
  ss.link <- get_ss_link(fit)
  traps <- get_trap(fit)[[session]]
  n.traps <- nrow(traps)
  
  # Perhaps need to add some validation for dists & mask
  if (is.null(dists)) {
    dists <- distances(traps, as.matrix(mask))
  }
  
  # A list to store each call density, will be combined into one data 
  density_list <- vector("list", length(call_id))
  
  for (i in seq_along(call_id)) {
    cid <- call_id[i]
    
    # Grab the capt & bincapt for this call ID
    sub_capt <- subset(capt, ID == cid)
    # Make sure to remove ID column when doing calculations
    sub_bin  <- subset(bincapt, ID == cid)[,-1] 
    ss.capt <- sub_capt$ss
    triggered_traps <- sub_capt$trap
    untriggered_traps <- setdiff(seq(n.traps), triggered_traps)
    
    # Initialize density matrix for this call
    mask.dens <- matrix(1, nrow = n.traps, ncol = nrow(mask))
    
    # Similar to how we calculate contribution from detection history alone
    for (j in seq_along(triggered_traps)) {
      mu.ss <- b0.ss - b1.ss*dists[triggered_traps[j], ]
      mask.dens[triggered_traps[j], ] <- dnorm(ss.capt[j], mu.ss, sigma.ss)
    }
    for (j in seq_along(untriggered_traps)) {
      mu.ss <- b0.ss - b1.ss*dists[untriggered_traps[j], ]
      mask.dens[untriggered_traps[j], ] <- pnorm(cutoff, mu.ss, sigma.ss)
    }
    
    mask.dens <- apply(mask.dens, 2, prod)
    
    density_list[[i]] <- data.frame(
      ID = cid, x = mask[, 1], y = mask[, 2], ss_density = mask.dens
    )
  }
  
  return(do.call(rbind, density_list))
}

#' Title
#'
#' @param fit 
#' @param call_id 
#' @param animal_id 
#' @param session 
#' @param mask 
#' @param dists 
#'
#' @return
#' @export
#'
#' @examples
capt_density <- function(fit, call_id, animal_id=NULL, session=1, 
                       mask=get_mask(fit)[[session]], dists=NULL) {
  # Note input validation for call_id, animal_id, session
  # is handled in `get_capt_by_id()`
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
  
  # Grab our estimated parameters, along with other model fitting parameters
  det.pars <- as.list(coef(fit, type="fitted"))
  if (fit$fit.types["ss"]){
    det.pars$cutoff <- get_ss.opts(fit)$cutoff
  }
  detfn <- get_detfn(fit)
  ss.link <- get_ss_link(fit)
  traps <- get_trap(fit)[[session]]
  
  # Perhaps need to add some validation for dists & mask
  if (is.null(dists)) {
    dists <- distances(traps, as.matrix(mask))
  }
  
  # A list to store each call density, will be combined into one data 
  density_list <- vector("list", length(call_id))
  
  for (i in seq_along(call_id)) {
    cid <- call_id[i]
    
    # Grab the capt & bincapt for this call ID
    sub_capt <- subset(capt, ID == cid)
    # Make sure to remove ID column when doing calculations
    sub_bin  <- as.numeric(subset(bincapt, ID == cid)[,-1]) 
    
    # Probability of detection at each mask point
    det.probs <- det_prob(detfn, det.pars, dists, ss.link)
    
    # Note that det_prob() returns E[SS] for ss models
    if(detfn == 'ss'){
      det.probs = 1 - pnorm(det.pars$cutoff, mean = det.probs, 
                            sd = det.pars$sigma.ss)
    }
    
    # Density due to capture history alone
    capt_dens <- apply((det.probs ^ sub_bin) * 
                         ((1 - det.probs) ^ (1 - sub_bin)), 2, prod)
    
    density_list[[i]] <- data.frame(
      ID = cid, x = mask[, 1], y = mask[, 2], capt_density = capt_dens
    )
  }
  
  return(do.call(rbind, density_list))
}

#' Title
#'
#' @param fit 
#' @param session 
#' @param mask 
#' @param call_id 
#'
#' @return
#' @export
#'
#' @examples
pois_density <- function(fit, call_id=NA, session=1, 
                         mask=get_mask(fit)[[session]], new_cov_data=NULL) {
  # Grab our estimated parameters, along with other model fitting parameters
  det.pars <- as.list(coef(fit, type="fitted"))
  if (fit$fit.types["ss"]){
    det.pars$cutoff <- get_ss.opts(fit)$cutoff
  }
  detfn <- get_detfn(fit)
  ss.link <- get_ss_link(fit)
  traps <- get_trap(fit)[[session]]

  a <- attr(mask, "area")
  
  ## Calculating density due to animal locations.
  p.det.total <- p.dot.defaultD(as.matrix(mask), traps, detfn, ss.link, det.pars, a)
  p.det <- p.dot.defaultD(as.matrix(mask), traps, detfn, ss.link, det.pars, a, 
                          returnNumerator = T)
  
  ## Divide by normalising constant; not conversion to square metres.
  if (fit$fit.ihd){
    if (is.null(new_cov_data)){
      D.mask.sess <- fit$D.mask[[session]]
    } else {
      D.mask.sess <- predict(fit, newdata=new_cov_data)
    }
    pois_dens <- (p.det * D.mask.sess) / (a * sum(p.det * D.mask.sess))
  } else {
    pois_dens <- p.det / p.det.total
  }

  return(data.frame(ID = call_id, 
                    x = mask[,1], 
                    y = mask[,2], 
                    pois_density = pois_dens))
}

#' Title
#'
#' @param call_id 
#' @param animal_id 
#' @param session 
#'
#' @return
#' @export
#'
#' @examples
create_locations_plot_subtitle <- function(call_id, animal_id, session) {
  if (is.null(animal_id)) {
    subtitle = paste0("session: ", session, ", call ID: ", call_id)
  } else {
    if (length(call_id) == 1) {
      subtitle = paste0("session: ", session, ", animal ID: ", animal_id,
                     ", call ID: ", call_id)
    } else {
      subtitle = paste0("session: ", session, ", animal ID: ", animal_id)
    }
  }
  
  return(subtitle)
}













