locations2 <- function(fit, 
                       id = "all", 
                       infotypes, 
                       session = 1, 
                       combine, 
                       xlim = NULL, 
                       ylim = NULL,
                       mask = get.mask(fit, session), 
                       nlevels = 10, 
                       cols = list(combined = 1, capt = 2, ss = 3, bearing = 4, dist = 5, toa = 6),
                       ltys = list(combined = "solid", capt = "solid",
                                   ss = "solid", bearing = "solid", dist = "solid", toa = "solid"),
                       trap.col = "red", circle.traps = TRUE,
                       show.labels = TRUE, plot.contours = TRUE,
                       plot.estlocs = FALSE, keep.estlocs = FALSE,
                       plot.arrows = "bearing" %in% fit$infotypes,
                       plot.circles = "dist" %in% fit$infotypes &
                         !("bearing" %in% fit$infotypes),
                       arrow.length = NULL, show.legend = FALSE,
                       show.axes = TRUE, add = FALSE, ask = TRUE) {
  # Check if animal_ID model
  animal.model = ('animal_ID' %in% colnames(fit$args$capt))
  
  ## Extracting the session's capture history.
  capt.all <- get.capt(fit, session)
  
  ## Setting id properly if "all" selected.
  if (id == "all"){
    ids <- 1:nrow(capt.all$bincapt) 
  }
  
  ## Sorting out limits.
  if (is.null(xlim)){
    xlim <- range(mask[, 1])
  }
  if (is.null(ylim)){
    ylim <- range(mask[, 2])
  }
  
  ## Setting up plotting area.
  if (!add){
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, asp = 1)
    box()
    if (show.axes){
      axis(1)
      axis(2)
    }
  }
  
  traps <- get.traps(fit, session)
  detfn <- get_detfn(fit)
  ss.link <- get_ss_link(fit)
  if (fit$fit.types["ss"]){
    fit$coefficients[["cutoff"]] <- fit$args$ss.opts$cutoff
  }
  dists <- distances(traps, mask)
  
  # Temporary ####
  det.probs <- det_prob(get_detfn(fit), fit$coefficients, dists, ss.link)
  p.det <- p.dot.defaultD(points = mask, traps = traps, detfn = detfn, pars=as.list(fit$coefficients), A=attr(mask, 'area'), ss.link)
  ################
  
  if(ask){
    ## Setting par(ask).
    ask.save <- par("ask")
    par(ask = TRUE)
    ## Making sure par is restored on function exit.
    on.exit(par(ask = ask.save))
  }
  
  plot.types <- c("combined", "capt", "ss", "bearing", "dist", "toa",  "mrds") %in% infotypes
  names(plot.types) <- c("combined", "capt", "ss", "bearing", "dist", "toa", "mrds")
  
  for (id in ids){
    densities <- calculate_densities(fit, id, infotypes)
    densities$location <- rep(1 / p.det, nrow(mask))
    
    if (plot.types["combined"]) {
      show.contour(mask = mask, 
                   dens = densities$combined, 
                   levels = NULL,
                   nlevels = nlevels, prob = FALSE, col = cols$capt,
                   lty = ltys$capt, show.labels = show.labels,
                   plot.contours = plot.contours, ask=ask)
    }
    # if (plot.types["capt"]) {
    #   if (!combine) {
    #     show.contour(mask = mask, 
    #                  dens = densities$capt * densities$location, 
    #                  levels = NULL,
    #                  nlevels = nlevels, prob = FALSE, col = cols$capt,
    #                  lty = ltys$capt, show.labels = show.labels,
    #                  plot.contours = plot.contours, ask=ask)
    #   }
    # }
    # if (plot.types["bearing"]) {
    #   if (!combine) {
    #     show.contour(mask = mask, 
    #                  dens = densities$bearing * densities$location, 
    #                  levels = NULL,
    #                  nlevels = nlevels, prob = FALSE, col = cols$capt,
    #                  lty = ltys$capt, show.labels = show.labels,
    #                  plot.contours = plot.contours, ask=ask)
    #   }
    # }
    # if (plot.types["dist"]) {
    #   if (!combine) {
    #     show.contour(mask = mask, 
    #                  dens = densities$dist * densities$location, 
    #                  levels = NULL,
    #                  nlevels = nlevels, prob = FALSE, col = cols$capt,
    #                  lty = ltys$capt, show.labels = show.labels,
    #                  plot.contours = plot.contours, ask=ask)
    #   }
    # }
    # if (plot.types["toa"]) {
    #   if (!combine) {
    #     show.contour(mask = mask, 
    #                  dens = densities$toa * densities$location, 
    #                  levels = NULL,
    #                  nlevels = nlevels, prob = FALSE, col = cols$capt,
    #                  lty = ltys$capt, show.labels = show.labels,
    #                  plot.contours = plot.contours, ask=ask)
    #   }
    # }
    
    capt <- capt.all$bincapt[1, ]
    ## Plotting traps, and circles around them.
    if (!add){
      points(traps, col = trap.col, pch = 4, lwd = 2)
      if (circle.traps){
        if (length(id) == 1){
          points(traps[capt == 1, , drop = FALSE], col = trap.col, cex = 2, lwd = 2)
        }
      }
    }
    ## Plotting arrows for estimated bearings.
    if (fit$fit.types["bearing"]){
      if (plot.arrows){
        if (fit$fit.types["dist"]){
          arrow.length <- capt.all$dist[id, capt == 1]
        }
        show.arrows(fit, id, session, arrow.length, trap.col)
      }
    }
    ## Plotting circles for estimated distances.
    if (fit$fit.types["dist"]){
      if (plot.circles & !plot.arrows){
        show.circles(fit, id, session, trap.col)
      }
    }
  }
}

#' Calculate the density of locations, capture history and covariates
#'
#' @param fit ACRE fit object
#' @param id
#' @param infotypes
#'
#' @return
#' @export
#'
#' @examples calculate_density(bd.fit, 2, infotypes = c("location", "capt", "dist", "bearing"))
calculate_densities <- function (fit, id, infotypes) {
  densities <- data.frame()

  # Calculate respective densities
  # Note that if an info-type was not used in the fitted model,
  # get_..._density(fit) will return 1
  location_density <- get_location_density(fit, id, infotypes)
  capt_density <- get_capt_density(fit, id, infotypes)
  dist_density <- get_dist_density(fit, id, infotypes)
  bearing_density <- get_bearing_density(fit, id, infotypes)
  ss_density <- get_ss_density(fit, id, infotypes)
  toa_density <- get_toa_density(fit, id, infotypes)

  # Note that if the density type is not included in 'infotypes' we are just
  # multiplying by 1
  combined_density <- location_density * capt_density * dist_density * 
    bearing_density * ss_density * toa_density
  
  densities <- list()
  densities$capt <- capt_density
  densities$location <- location_density
  densities$dist <- dist_density
  densities$bearing <- bearing_density
  densities$ss <- ss_density
  densities$toa <- toa_density
  densities$combined <- combined_density
  # densities <- data.frame(
  #   location = location_density,
  #   capt = capt_density,
  #   dist = dist_density,
  #   bearing = bearing_density,
  #   ss = ss_density,
  #   toa = toa_density,
  #   combined = combined_density
  # )
  
  return(densities)
}

get_location_density <- function(fit, id, infotypes) {
  if ("location" %in% infotypes) {
    # return(fit$output.tmb$report$location_density)
    return(1)
  } else {
    return(1)
  }
}

get_capt_density <- function(fit, id, infotypes) {
  if ("capt" %in% infotypes) {
    # return(fit$output.tmb$report$capt_density)
    return(fit$output.tmb$report$density_array[id, 1, ])
  } else {
    return(1)
  }
}

get_dist_density <- function(fit, id, infotypes) {
  if ("dist" %in% infotypes) {
    # return(fit$output.tmb$report$dist_density)
    return(fit$output.tmb$report$density_array[id, 2, ])
  } else {
    return(1)
  }
}

get_bearing_density <- function(fit, id, infotypes) {
  if ("bearing" %in% infotypes) {
    # return(fit$output.tmb$report$bearing_density)
    return(fit$output.tmb$report$density_array[id, 3, ])
  } else {
    return(1)
  }
}

get_ss_density <- function(fit, id, infotypes) {
  if ("ss" %in% infotypes) {
    # return(fit$output.tmb$report$ss_density)
    return(fit$output.tmb$report$density_array[id, 4, ])
  } else {
    return(1)
  }
}

get_toa_density <- function(fit, id, infotypes) {
  if ("toa" %in% infotypes) {
    # return(fit$output.tmb$report$toa_density)
    return(fit$output.tmb$report$density_array[id, 5, ])
  } else {
    return(1)
  }
}







