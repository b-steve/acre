#' Plotting estimated locations
#'
#' Plots estimated densities of animal locations, which are latent
#' variables in SECR models.
#'
#' @param fit A fitted model from \link{fit.ascr}.
#' @param id A numeric vector with row numbers from
#'     \code{fit$args$capt}, indicating which individuals' locations
#'     are to be plotted. Alternatively, the character string
#'     \code{"all"}, indicating all animals within the selected
#'     session.
#' @param session The session with the detector array and invidual(s)
#'     to be plotted (for multi-session models only).
#' @param infotypes A character vector indicating the type(s) of
#'     information to be used when plotting the estimated density of
#'     location.  Elements can be a subset of \code{"capt"},
#'     \code{"bearing"}, \code{"dist"}, \code{"ss"}, \code{"toa"},
#'     \code{"combined"}, and \code{"all"}, where \code{"capt"} shows
#'     estimated location only using detection locations,
#'     \code{"combined"} combines all information types together, and
#'     \code{"all"} plots all possible contour types. When signal
#'     strength information is used in the model fit, \code{"capt"}
#'     and \code{"ss"} are equivalent as the signal strength
#'     information is built into the detection function. By default,
#'     only the most informative contour is plotted, i.e.,
#'     \code{"capt"} if the model was fitted with no additional
#'     information, and \code{"combined"} otherwise.
#' @param combine Logical, if \code{TRUE} then the information types
#'     specified in \code{infotypes} are combined into a single
#'     contour. If \code{FALSE} then separate contours are plotted for
#'     each information type.
#' @param xlim A numeric vector of length 2, giving the x coordinate
#'     range.
#' @param ylim A numeric vector of length 2, giving the y coordinate
#'     range.
#' @param mask A matrix with two columns. Each row provides Cartesian
#'     coordinates for the location of a mask point. The function
#'     \link[ascr]{create.mask} will return a suitable object. The
#'     mask used to fit the model \code{fit} will be used by default;
#'     this argument is usually used when estimated location contours
#'     need to be plotted to a higher resolution than this.
#' @param newdata If a new mask is provided in the \code{mask}
#'     argument, and the model fits an inhomogeneous density surface
#'     via mask-level covariates, then this argument must provide the
#'     covariates for the new mask's point locations.
#' @param levels A numeric vector giving the values to be associated
#'     with the plotted contours.
#' @param nlevels The number of contour levels desired. Ignored if
#'     \code{levels} is provided.
#' @param density Logical, if \code{TRUE}, the labels on contours (and
#'     the levels specified by \code{levels}) refer to the density of
#'     the estimated distribution of the individual's location. If
#'     \code{FALSE}, the labels on contours (and the levels specified
#'     by \code{levels}) refer to the probability of the individual
#'     being located within the associated contour under the estimated
#'     distribution of the individual's location.
#' @param cols A list with named components corresponding to each
#'     contour type (i.e., a subset of \code{"capt"},
#'     \code{"bearing"}, \code{"dist"}, \code{"toa"}, and
#'     \code{"combined"}). Each component provides the colour of the
#'     associated contour type (e.g., using a character string such as
#'     \code{"red"}, or a call to the function
#'     \link[grDevices]{rgb}). By default, if only one contour is to
#'     be plotted, it will be plotted in black. Alternatively, a
#'     vector with a single element, specifying the colour for all
#'     contours.
#' @param ltys The line type of the contours, with the same required
#'     syntax as \code{cols}; see \link{par}.
#' @param trap.col The colour of the points representing detector
#'     locations.
#' @param circle.traps Logical, if \code{TRUE} circles are plotted
#'     around traps that made a detection of the individual in
#'     question.
#' @param show.labels Logical, if \code{TRUE}, contours are labelled
#'     with the appropriate probability density (if \code{density} is
#'     \code{TRUE}), or the corresponding probability of the
#'     individual being within the associated contour, under the
#'     estimated density (if \code{density} is \code{FALSE}).
#' @param plot.contours Logical, if \code{TRUE}, contours are
#'     plotted. Note that, if \code{FALSE}, nothing corresponding to
#'     the density of the individuals' locations is plotted unless
#'     \code{plot.estlocs} is \code{TRUE}.
#' @param plot.estlocs Logical, if \code{TRUE}, dots are plotted at
#'     the mode of the combined densities. If a density has more than
#'     a single mode (and the modes have the same density value) then
#'     a dot will only be plotted at one of them.
#' @param keep.estlocs Logical, if \code{TRUE}, the locations of the
#'     estimated locations are invisibly returned.
#' @param plot.arrows Logical, if \code{TRUE}, arrows indicating the
#'     estimated bearing to the individual are plotted from detectors
#'     at which detections were made.
#' @param plot.circles Logical, if \code{TRUE}, circles indicating the
#'     estimated distance to the individual are plotted around
#'     detectors at which detections were made.
#' @param arrow.length Numeric, providing the length of the arrows
#'     (only used if \code{plot.arrows} is \code{TRUE}).
#' @param show.legend Logical, if \code{TRUE}, a legend will be added
#'     to the plot.
#' @param show.axes Logical, if \code{TRUE}, axes will be added to the
#'     plot.
#' @param add Logical, if \code{TRUE}, contours will be added to an
#'     existing plot.
#' @param ask Logical. If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn.
#'
#' @return If \code{keep.estlocs} is \code{TRUE}, then a list
#'     containing a matrix of estimated locations is invisibly
#'     returned. See examples.
#' 
#' @examples
#' locations(example.data$fits$simple.hn, 1)
#' locations(example.data$fits$simple.hn, 1, levels = c(0.50, 0.90, 0.95))
#' ## Saving estimated locations.
#' estlocs <- locations(example.data$fits$simple.hn, keep.estlocs = TRUE)
#' estlocs
#' \dontrun{
#' fine.mask <- create.mask(example.data$traps, 20, spacing = 0.2)
#' locations(example.data$fits$bearing.hn, 1, infotypes = "all", mask = fine.mask)
#' }
#'
#' @export
locations <- function(fit, call_id = "all", animal_id=NULL, session = 1, infotypes = NULL,
                      combine = FALSE, xlim = NULL, ylim = NULL,
                      mask = get_mask(fit)[[session]], newdata = NULL,
                      levels = NULL, nlevels = 10, density = FALSE,
                      cols = list(combined = 1, capt = 2, ss = 3, bearing = 4, dist = 5, toa = 6),
                      ltys = list(combined = "solid", capt = "solid",
                                  ss = "solid", bearing = "solid", dist = "solid", toa = "solid"),
                      trap.col = "red", circle.traps = TRUE,
                      show.labels = TRUE, plot.contours = TRUE,
                      plot.estlocs = FALSE, keep.estlocs = FALSE,
                      plot.arrows = "bearing" %in% fit$infotypes,
                      plot.ss = "ss" %in% fit$infotypes,
                      plot.circles = "dist" %in% fit$infotypes &
                          !("bearing" %in% fit$infotypes),
                      arrow.length = NULL, show.legend = FALSE,
                      show.axes = TRUE, add = FALSE, ask = TRUE, 
                      joint.density = FALSE) {
  
    ## Type checking
    stopifnot(
      # Check that 'session' is numeric and has length 1
      is.numeric(session), 
      length(session) == 1,
      
      # Check that 'animal_id' is either NULL or a single numeric value
      is.null(animal_id) || (is.numeric(animal_id) && length(animal_id) == 1),
      
      # Check that 'call_id' is either "all", a single numeric value, or a numeric vector
      (is.character(call_id) && call_id == "all") || 
        (is.numeric(call_id) && all(!is.na(call_id)))
    )
  
    ## Error if session argument is too large.
    if (session > fit$n.sessions){
        if (fit$n.sessions == 1){
            stop(paste("Argument 'session' too large; there was only ", fit$n.sessions, " session.", sep = ""))
        } else {
            stop(paste("Argument 'session' too large; there were only ", fit$n.sessions, " sessions.", sep = ""))
        }
    }
  
    ## Error for locations() with a directional model.
    if (!is.null(fit$args$ss.opts$directional)){
        if (fit$args$ss.opts$directional){
            stop("The locations() function has not yet been implemented for directional model fits.")
        }
    }
    if (!is.null(fit$args$ss.opts$het.source)){
        if (fit$args$ss.opts$het.source){
            stop("The locations() function has not yet been implemented for heterogeneous source strength model fits.")
        }
    }
    if (!is.null(fit$first.calls)){
        if (fit$first.calls){
            stop("The locations() function has not yet been implemented for first-call models.")
        }
    }
    ## Error if combine specified without infotypes.
    if (missing(infotypes) && combine){
        stop("Argument `combine' is only useful if `infotypes' is provided.")
    }
    ## Error if new mask covariates are not provided when they need to be.
    is.new.mask <- !missing(mask)
    if (fit$fit.ihd & is.new.mask & is.null(newdata)){
        stop("Covariate values for the mask object must be provided via the `newdata' argument.")
    }
    
    ## Extracting the session's capture history.
    # capt.all <- get_capt(fit, session)
    # TODO: Multi session
    # capt.all <- get_capt_for_plot(fit$args)
    # capt.all <- capt.all[capt.all$session == session, ]
    
    capt.all <- fit$arg_input$captures
    capt.all <- capt.all[capt.all$session == session, ]
    
    animal.model <- is_animal_model(fit)
    
    traps <- get_trap(fit)[[session]]
    
    if(animal.model && is.null(animal_id)) {
      stop("'animal_id' must be provided for animal ID models")
    }
    
    if(!is.null(animal_id) && !animal.model) {
      stop("'animal_id' argument should only be provided for animal models. Use 'call_id' argument instead.")
    }
    
    if (!is.null(animal_id)) {
      if (!(animal_id) %in% capt.all$animal_ID) {
        stop(paste("Could not find 'animal_id'", animal_id))
      }
      
      capt.all <- subset(capt.all, capt.all$animal_ID == animal_id)
    }
    
    ## Setting id properly if "all" selected.
    if (any(call_id == "all")) {
      call_id <- unique(capt.all$ID)
    }
    
    if (joint.density) {
      message("Calculating animal densities: ")
      progress_bar = txtProgressBar(min = 0, max = length(call_id), initial = 0, style=3) 
    }

    ## Saving estimated locations.
    if (keep.estlocs){
        estlocs <- matrix(0, nrow = length(call_id), ncol = 2)
        j <- 1
    }
    
    ## Sorting out limits.
    
    stopifnot(all(is.numeric(xlim) || is.null(xlim), 
                  is.numeric(ylim) || is.null(ylim)))
    
    if (is.null(xlim)){
        xlim <- range(mask[,1])
    }
    if (is.null(ylim)){
        ylim <- range(mask[,2])
    }
    
    
    # Setup plotting canvas
    base.plot <- base_plot(capt.all, as.data.frame(mask), xlim, ylim) + scale_shape_identity()
    
    
    ## Ignoring 'nlevels' if 'levels' is provided.
    if (!is.null(levels)){
        nlevels <- length(levels)
    }
    
    ## Logical value for whether or not any additional information was
    ## used in model fit.
    any.infotypes <- length(fit$infotypes) > 0
    ## Setting default infotypes.
    if (is.null(infotypes)){
        if ("mrds" %in% fit$infotypes){
            infotypes <- "mrds"
        } else {
            if (any.infotypes){
                infotypes <- "combined"
            } else {
                infotypes <- "capt"
            }
        }
    }
    
    warn.estlocs <- FALSE
    if (any.infotypes){
        if (!("combined" %in% infotypes) & (plot.estlocs | keep.estlocs)){
            warn.estlocs <- TRUE
        }
    } else {
        if (!("capt" %in% infotypes) & (plot.estlocs | keep.estlocs)){
            warn.estlocs <- TRUE
        }
    }
    if (warn.estlocs){
        warning("Estimated locations only available for 'combined' information type. Ignoring 'plot.estlocs' and 'keep.estlocs'.")
    }
    ## Removing "all" or "combined" from infotypes if "mrds" is specified.
    if ("mrds" %in% infotypes){
        infotypes <- infotypes[infotypes != "all"]
        infotypes <- infotypes[infotypes != "combined"]
    }
    ## Error if "combined" is used when there is no additional information.
    if (("combined" %in% infotypes) & !any.infotypes){
        stop("No additional information used in model 'fit', so a \"combined\" contour cannot be plotted.")
    }
    ## Working out which contours to plot.
    if ("all" %in% infotypes){
        # TODO: is this supposed to be like this?
        infotypes <- c(fit$infotypes, "capt", "combined"[any.infotypes])
    }
    
    ## If "ss" is an infotype, set to "capt". OR NOT. <<- Wait what the hell is going on in this comment. <<- lmao
    ##infotypes[infotypes == "ss"] <- "capt"
    
    infotypes <- unique(infotypes)
    
    ## Setting colour to "black" if there is only one contour to be plotted.
    if (missing(cols)) {
        if (length(infotypes) == 1){
            cols <- "black"
        }
    }
    if (missing(ltys)){
        if (length(infotypes) == 1){
            ltys <- "solid"
        }
    }
    if (is.list(cols)){
        if (!all(infotypes %in% names(cols))){
            stop("Provide a colour for each contour to be plotted.")
        }
    }
    if (is.list(ltys)){
        if (!all(infotypes %in% names(ltys))){
            stop("Provide a line type for each contour to be plotted.")
        }
    }
    if (length(cols) == 1){
        cols.save <- cols
        cols <- vector(mode = "list", length = length(infotypes))
        names(cols) <- infotypes
        cols[infotypes] <- cols.save
    }
    if (length(ltys) == 1){
        ltys.save <- ltys
        ltys <- vector(mode = "list", length = length(infotypes))
        names(ltys) <- infotypes
        ltys[infotypes] <- ltys.save
        if (combine){
            ltys["combined"] <- ltys.save
        }
    }
    plot.types <- c("combined", "capt", "ss", "bearing", "dist", "toa",  "mrds") %in% infotypes
    names(plot.types) <- c("combined", "capt", "ss", "bearing", "dist", "toa", "mrds")
    if (combine){
        plot.types["combined"] <- TRUE
    }
    ## Setting all to TRUE if combined is used.
    ## Some error catching.
    for (i in c("bearing", "dist", "toa")){
        if (plot.types[i] & !fit$fit.types[i]){
            msg <- paste("Contours for information type '", i, "' cannot be plotted as this information was not used in the model 'fit'", sep = "")
            warning(msg)
            plot.types[i] <- FALSE
        }
    }
    det.pars <- as.list(coef(fit, type="fitted"))
    if (fit$fit.types["ss"]){
      det.pars$cutoff <- get_ss.opts(fit)$cutoff
    }
    detfn <- get_detfn(fit)
    ss.link <- get_ss_link(fit)
    dists <- distances(traps, mask)
    a <- attr(mask, "area")
    
    ## Calculating density due to animal locations.
    p.det <- p.dot.defaultD(mask, traps, detfn, ss.link, det.pars, a)
    
    ## Divide by normalising constant; not conversion to square metres.
    if (fit$fit.ihd){
        if (is.new.mask){
            D.mask.sess <- predict(fit, newdata)
        } else {
            D.mask.sess <- fit$D.mask[[session]]
        }
        f.x <- p.det*D.mask.sess/(a*sum(p.det*D.mask.sess))
    } else {
        f.x <- p.det/(a*sum(p.det))
    }
    
    n_mask_points <- nrow(mask)
    densities.list <- list(
      f.capt = rep(1, n_mask_points),
      f.ss = rep(1, n_mask_points),
      f.bearing = rep(1, n_mask_points),
      f.dist = rep(1, n_mask_points),
      f.toa = rep(1, n_mask_points),
      f.combined = rep(1, n_mask_points)
    )
    
    for (i in call_id){
        # Sort out plot title
        if (animal.model) {
          if (joint.density) {
            sub_title = paste0("session: ", session, ", animal ID: ", i)
          } else {
            sub_title = paste0("session: ", session, ", animal ID: ", animal_id,
                               ", call ID: ", i)
          }
        } else {
          sub_title = paste0("session: ", session, ", call ID: ", i)
        }
        
        base.plot <- base.plot + labs(subtitle = sub_title)
        
        if (plot.types["combined"]){
            if ((!combine) | (combine & plot.types["capt"])){
                f.combined <- f.x
                densities.list$f.combined <- densities.list$f.combined * f.x
                
            } else {
                f.combined <- 0*f.x + 1
            }
        }
        
        capt <- get_capt_by_id(fit, i, animal_id, session)
        bincapt <- get_capt_by_id(fit, i, animal_id, session, return_bincapt=T)
        
        ## Contour due to capture history.
        if (plot.types["capt"] | plot.types["combined"] | plot.types["ss"]){
            det.probs <- det_prob(detfn, det.pars, dists, ss.link)
            
            if(detfn == 'ss'){
              det.probs = 1 - pnorm(det.pars$cutoff, mean = det.probs, sd = det.pars$sigma.ss)
            }

            f.capt <- apply((det.probs ^ bincapt) * ((1 - det.probs) ^ (1 - bincapt)), 2, prod)
            
            # Update joint density
            densities.list$f.capt <- densities.list$f.capt * f.capt
            
            if (plot.types["capt"]){
                if (!combine & !joint.density){
                  contours <- calculate.contour(mask = mask, dens = f.x*f.capt, levels = levels,
                                 nlevels = nlevels, prob = !density)
                  
                  capt_normalised_density <- f.x*f.capt / (attr(mask, "a")*10000*sum(f.x*f.capt)) 
                  
                  capt_contour <- ggplot2::geom_contour(
                    data = data.frame(x=mask[,1], y=mask[,2], z=capt_normalised_density),
                    aes(x = x, y = y, z = z),
                    breaks=contours$levels,
                    color = cols$capt,
                    linetype= ltys$capt
                    )
                  
                  base.plot <- base.plot + capt_contour
                }
            }
            if (fit$fit.types["ss"]){
                f.ss.capt <- ss.density(fit, session, mask, dists, call_id = i, animal_id = animal_id)
                f.ss <- f.ss.capt/f.capt
                
                # Update joint density
                densities.list$f.ss <- densities.list$f.ss * f.ss
                
                ## Such a hack, but this keeps f.combined correct,
                ## below.
                f.capt <- f.ss.capt
                ## Ughhhh sorry about this one.
                f.ss[f.ss == Inf] <- 0
                if (plot.types["ss"]){
                    if (!combine & !joint.density){
                      contours <- calculate.contour(mask = mask, dens = f.x*f.ss, levels = levels,
                                     nlevels = nlevels, prob = !density)
                      
                      ss_normalised_density <- t(f.x*f.ss / (attr(mask, "a")*10000*sum(f.x*f.ss)))
                      
                      ss_contour <- ggplot2::geom_contour(
                        data = data.frame(x=mask[,1], y=mask[,2], z=t(ss_normalised_density)),
                        aes(x = x, y = y, z = z),
                        breaks=contours$levels,               
                        color = cols$ss,
                        linetype= ltys$ss
                      )
                      
                      base.plot <- base.plot + ss_contour
                    }
                }
            } else {
                f.ss <- f.capt*0 + 1
            }
            ## This shit is such a mess, sorry if you have to work out
            ## how this works later.
            if (plot.types["combined"] & !combine){
                f.combined <- f.combined*f.capt
                densities.list$f.combined <- densities.list$f.combined * f.capt
            } else if (plot.types["combined"] & combine){
                f.true.capt <- f.capt/f.ss
                f.true.capt[f.ss == 0] <- 0
                if (plot.types["capt"] & plot.types["ss"]){
                    f.combined <- f.combined*f.capt
                    densities.list$f.combined <- densities.list$f.combined * f.capt
                } else if (plot.types["capt"] & !plot.types["ss"]){
                    f.combined <- f.combined*f.true.capt
                    densities.list$f.combined <- densities.list$f.combined * f.true.capt
                } else if (!plot.types["capt"] & plot.types["ss"]){
                    f.combined <- f.combined*f.ss
                    densities.list$f.combined <- densities.list$f.combined * f.ss
                }
            }
        }
        ## Contour due to estimated bearings.
        if (plot.types["bearing"] | plot.types["combined"] & fit$fit.types["bearing"]){
            f.bearing <- bearing.density(fit, session, mask, call_id = i, animal_id = animal_id)
            
            # Update joint density
            densities.list$f.bearing <- densities.list$f.bearing * f.bearing
            
            if (plot.types["bearing"]) {
                if (!combine & !joint.density) {
                  contours <- calculate.contour(mask = mask, dens = f.x*f.bearing, levels = levels,
                                 nlevels = nlevels, prob = !density)
                  
                  bearing_normalised_density <- t(f.x*f.bearing / (attr(mask, "a")*10000*sum(f.x*f.bearing)))
                  
                  bearing_contour <- ggplot2::geom_contour(
                    data = data.frame(x=mask[,1], y=mask[,2], z=bearing_normalised_density),
                    aes(x = x, y = y, z = z),
                    breaks=contours$levels,                 
                    color = cols$bearing,
                    linetype= ltys$bearing
                  )

                  base.plot <- base.plot + bearing_contour
                }
            }
            if (plot.types["combined"]){
                if ((!combine) | (combine & plot.types["bearing"])){
                    f.combined <- f.combined*f.bearing
                    densities.list$f.combined <- densities.list$f.combined * f.bearing
                }
            }
        }
        ## Contour due to estimated distances.
        if (plot.types["dist"] | plot.types["combined"] & fit$fit.types["dist"]){
            f.dist <- dist.density(fit, session, mask, dists, call_id = i, animal_id = animal_id)
            
            # Update joint density
            densities.list$f.dist <- densities.list$f.dist * f.dist
            
            if (plot.types["dist"]){
                if (!combine & !joint.density){
                  contours <- calculate.contour(mask = mask, dens = f.x*f.dist, levels = levels,
                                 nlevels = nlevels, prob = !density)
                  
                  dist_normalised_density <- t(f.x*f.dist / (attr(mask, "a")*10000*sum(f.x*f.dist)))
                  
                  dist_contour <- ggplot2::geom_contour(
                    data = data.frame(x=mask[,1], y=mask[,2], z=dist_normalised_density),
                    aes(x = x, y = y, z = z),
                    breaks=contours$levels,  
                    color = cols$dist,
                    linetype= ltys$dist
                  )
                  
                  base.plot <- base.plot + dist_contour
                }
            }
            if (plot.types["combined"]){
                if ((!combine) | (combine & plot.types["dist"])){
                    f.combined <- f.combined*f.dist
                    densities.list$f.combined <- densities.list$f.combined * f.dist
                }
            }
        }
        ## Contour due to measured times of arrival.
        # Note toa data not useful if only 1 trap is triggered
        if ((plot.types["toa"] | plot.types["combined"]) &
            fit$fit.types["toa"] & sum(bincapt) > 1){
            f.toa <- toa.density(fit, session, mask, dists, call_id = i, animal_id = animal_id)
            
            # Update joint density
            densities.list$f.toa <- densities.list$f.toa * f.toa
            
            if (plot.types["toa"]){
                if (!combine & !joint.density){
                  contours <- calculate.contour(mask = mask, dens = f.x*f.toa, levels = levels,
                                 nlevels = nlevels, prob = !density)
                  
                  toa_normalised_density <- f.x*f.toa / (attr(mask, "a")*10000*sum(f.x*f.toa))
                  
                  toa_contour <- ggplot2::geom_contour(
                    data = data.frame(x=mask[,1], y=mask[,2], z=toa_normalised_density),
                    aes(x = x, y = y, z = z),
                    breaks=contours$levels,
                    color = cols$toa,
                    linetype= ltys$toa
                  )
                  
                  base.plot <- base.plot + toa_contour
                }
            }
            if (plot.types["combined"]){
                  if ((!combine) | (combine & plot.types["toa"])){
                      f.combined <- f.combined*f.toa
                      densities.list$f.combined <- densities.list$f.combined * f.toa
                  }
              }
        }
        
        # Update joint density
        # densities.list$f.combined <- densities.list$f.combined * f.combined
        
        ## Combined contour.
        if (plot.types["combined"] & !joint.density){
            contours <- calculate.contour(mask = mask, dens = f.combined, levels = levels,
                         nlevels = nlevels, prob = !density)
            
            combined_normalised_density <- f.combined / (attr(mask, "a")*10000*sum(f.combined)) 

            base.plot <- plot.contour(
              base_plot = base.plot,
              x = mask[,1], y = mask[,2], 
              density = combined_normalised_density,
              levels = contours$levels,
              ltys = ltys$combined,
              color = NA
            )
        }
        if (plot.types["mrds"]){
            loc <- capt.all$mrds[call_id, , drop = FALSE]
            points(loc, pch = 16, col = "black")
        }
        if (plot.estlocs | keep.estlocs){
            if ((any.infotypes & plot.types["combined"]) | !any.infotypes){
                f.estlocs <- if (any.infotypes) f.combined else f.capt*f.x
                mode.points <- which(f.estlocs == max(f.estlocs))[1]
                if (plot.estlocs){
                    points(mask[mode.points, 1], mask[mode.points, 2],
                           pch = 16, col = "black")
                }
                if (keep.estlocs){
                    estlocs[j, ] <- c(mask[mode.points, 1], mask[mode.points, 2])
                    j <- j + 1
                }
            }
        }

        
        
        ## Plotting traps, and circles around them.
        if (!add){
          if (plot.ss & !joint.density) {
            plot_traps <- data.frame(traps[which(bincapt == 0), , drop = FALSE])
          } else {
            plot_traps <- data.frame(traps)
          }
          # Plot traps
          base.plot = base.plot + 
            geom_point(data = plot_traps, mapping = aes(x = x, y = y, shape = 4), 
                       size = 3, col="red", stroke=1)
            
          
          if ((circle.traps & !plot.ss) | (circle.traps & joint.density)){
            base.plot = base.plot + geom_point(
              data = subset(traps, bincapt == 1),
              aes(x = x, y = y),
              shape = 21,
              color = "red",
              size = 4.5,
              stroke = 1
            )
          }
          
          # Plotting arrows for estimated bearings.
          if (fit$fit.types["bearing"]) {
            if (plot.arrows) {
              base.plot = base.plot + estimated_bearing_arrow_plot(fit = fit,
                                                                   session = session,
                                                                   call_id = i,
                                                                   animal_id = animal_id)
            }
          }
          
          
          # Plotting distance circles
          if (fit$fit.types["dist"]) {
            if (plot.circles & !plot.arrows){
              base.plot = base.plot + estimated_distance_plot(fit = fit,
                                                              session = session,
                                                              call_id = i,
                                                              animal_id = animal_id)
            }
          }
          
          # Plotting ss indicator circles
          if (fit$fit.types["ss"]) {
            if (plot.ss & !joint.density){
              activated_traps <- traps[which(bincapt == 1), , drop = FALSE]
              base.plot = base.plot + geom_point(data = activated_traps, 
                                                 mapping = aes(x = activated_traps[, 1],
                                                               y = activated_traps[, 2],
                                                               colour = capt$ss),
                                                 size = 5)
            }
          }
        }
        
        
        if (!joint.density) {
          plot(base.plot)
        }
        
        if(!add & !joint.density) {
          # Reset the plot
          base.plot = base_plot(capt.all, as.data.frame(mask), xlim, ylim) + 
            scale_shape_identity()
        }
        
        if (ask & !joint.density) {
          # Make sure to only ask if we are plotting more than one plot, 
          # and it is not the last plot in the list
          if (length(call_id) > 1 && i != call_id[[length(call_id)]]) {
            prompt_user_for_next_plot()
          }
        }
        
        if (joint.density) {
          setTxtProgressBar(progress_bar, which(call_id == i))
        }
    }
    
    if (joint.density) {
      contours <- calculate.contour(mask = mask, dens = densities.list$f.combined, levels = levels,
                                    nlevels = nlevels, prob = !density)
      
      combined_normalised_density <- densities.list$f.combined / (attr(mask, "a")*10000*sum(densities.list$f.combined)) 
      # geom_contour()
      combined_contour <- ggplot2::geom_contour_filled(
        data = data.frame(x=mask[,1], y=mask[,2], z=combined_normalised_density),
        aes(x = x, y = y, z = z),
        breaks=contours$levels,
        # color = cols$combined,
        color = NA,
        linetype= ltys$combined
      )
      
      
      ### TESTING CUSTOM FILL LEVELS ###
      # Now define a manual fill scale, making the first bin transparent
      # You need one color for each bin. If you have length(contours$levels) - 1 bins,
      # you must supply that many values in scale_fill_manual().
      num_bins <- length(contours$levels) - 1
      

      # Create an interpolating function from the max number of colors 
      # that PuRd supports (which is 9).
      palette_fun <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
      
      # base_colors <- palette_fun(num_bins)
      base_colors <- viridisLite::viridis(num_bins)
      alpha_values <- seq(0, 1, length.out = num_bins)
      faded_colors <- unlist(
        Map(function(col, a) scales::alpha(col, a), base_colors, alpha_values)
      )
      
      # Now generate as many colors as needed:
      my_colors <- palette_fun(num_bins)
      
      # Prepend 'transparent' (or NA) for the lowest bin
      fill_values <- c("transparent", my_colors)
      
      ######
      
      base.plot <- base.plot + combined_contour + ggplot2::scale_fill_manual(values = unname(faded_colors))
      # base.plot <- base.plot + combined_contour +
      #   ggplot2::scale_fill_manual(values = c("transparent", base_colors)) + 
      # ggplot2::theme(
      #     panel.background = ggplot2::element_rect(fill = viridisLite::viridis(1), color = NA),
      #     panel.grid.major = ggplot2::element_blank(),
      #     panel.grid.minor = ggplot2::element_blank()
      #     )
          
      # base.plot <- base.plot + combined_contour
      
      close(progress_bar)
      
      plot(base.plot)
    }
    
    ## Making legend.
    if (show.legend){
      legend.labels <- infotypes
      legend.cols <- c(cols[infotypes], recursive = TRUE)
      legend.ltys <- c(ltys[infotypes], recursive = TRUE)
      legend("topright", legend = infotypes, lty = legend.ltys, col = legend.cols, bg = "white")
    }
    out <- list()
    if (keep.estlocs){
      out$estlocs <- estlocs
    }
    
    else {
      out <- TRUE
    }
    
    invisible(out)
}

plot.contour <- function(base_plot, x, y, density, levels, ltys, color=NA) {
  combined_contour <- ggplot2::geom_contour_filled(
    data = data.frame(x=x, y=y, z=density),
    aes(x = x, y = y, z = z),
    breaks=levels,
    # color = cols$combined,
    color = NA,
    linetype= ltys
  )
  
  num_bins <- length(levels) - 1
  
  
  # Create an interpolating function from the max number of colors 
  # that PuRd supports (which is 9).
  # palette_fun <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
  
  # base_colors <- palette_fun(num_bins)
  base_colors <- viridisLite::viridis(num_bins)
  alpha_values <- seq(0, 1, length.out = num_bins)
  faded_colors <- unlist(
    Map(function(col, a) scales::alpha(col, a), base_colors, alpha_values)
  )
  
  # Now generate as many colors as needed:
  # my_colors <- palette_fun(num_bins)
  
  # Prepend 'transparent' (or NA) for the lowest bin
  # fill_values <- c("transparent", my_colors)
  
  ######
  
  return(base_plot + 
           combined_contour + 
           ggplot2::scale_fill_manual(values = unname(faded_colors)))
}

## Function calculates contour levels based by probability rather than density
#'
#' @export
calculate.contour <- function(mask, dens, nlevels, levels, prob){
    ## Divide densities by normalising constant before plotting.
    a <- attr(mask, "a")*10000
    ## Note conversion of area to square metres.
    dens <- dens/(a*sum(dens))
    unique.x <- sort(unique(mask[, 1]))
    unique.y <- sort(unique(mask[, 2]))
    z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
    n.mask <- nrow(mask)
    for (i in 1:n.mask){
      x <- mask[i, 1]
      y <- mask[i, 2]
      index.x <- which(x == unique.x)
      index.y <- which(y == unique.y)
      z[index.x, index.y] <- dens[i]
    }
    ## Sorting out levels.
    if (is.null(levels)){
      levels <- pretty(range(z, finite = TRUE), nlevels)
    } else {
      if (prob){
        z.sort <- sort(z, decreasing = TRUE)
        probs.sort <- cumsum(z.sort)/sum(z.sort)
        prob.levels <- levels
        levels <- numeric(nlevels)
        for (i in 1:nlevels){
          levels[i] <- z.sort[which(abs(probs.sort - prob.levels[i]) ==
                                      min(abs(probs.sort - prob.levels[i])))[1]]
        }
      }
    }
    if (prob){
      labels <- character(nlevels)
      for (i in 1:nlevels){
        labels[i] <- format(round(sum(z[z > levels[i]], na.rm = TRUE)/
                                    sum(z, na.rm = TRUE), 2), nsmall = 2)
      }
    } else {
      labels <- NULL
    }
    
    return(list(x = unique.x, y = unique.y, z = z, levels = levels, 
                labels = labels))
}


## Calculating density due to estimated bearings.
bearing.density <- function(fit, session, mask, call_id, animal_id=NULL){
    if (is_animal_model(fit) & is.null(animal_id)) {
      stop("'animal_id' argument must be provided for animal id models")
    }
    
    capt <- get_capt_by_id(fit, call_id, animal_id, session)
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
    bearing.capt <- capt$bearing
    
    # Initialize density matrix for this call
    mask.dens <- matrix(0, nrow = sum(bincapt), ncol = nrow(mask))
    
    # Grab our estimated kappa parameter
    kappa <- coef(fit, type="fitted")[["kappa"]]
    
    # Calculate bearings for each triggered trap
    mask.bearings <- bearings(get_trap(fit)[[session]][bincapt == 1, , drop = FALSE], as.matrix(mask))
    
    # Calculate the density for each of our detections
    for (j in 1:sum(bincapt)) {
      mask.dens[j, ] <- CircStats::dvm(bearing.capt[j], mu = mask.bearings[j, ], kappa = kappa)
    }

    return(apply(mask.dens, 2, prod))
}

## Calculating density due to estimated distances.
dist.density <- function(fit, session, mask, dists, call_id, animal_id=NULL){
    if (is_animal_model(fit) & is.null(animal_id)) {
      stop("'animal_id' argument must be provided for animal id models")
    }
    
    capt <- get_capt_by_id(fit, call_id, animal_id, session)
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)

    dist.capt <- capt$dist
    
    # Initialize density matrix for this call
    mask.dens <- matrix(0, nrow = sum(bincapt), ncol = nrow(mask))
    
    # Grab our estimated alpha parameter
    alpha <- coef(fit, type="fitted")[["alpha"]]
    betas <- alpha/dists
    
    # Calculate bearings for each triggered trap
    for (j in 1:sum(bincapt)) {
      mask.dens[j, ] <- dgamma(dist.capt[j], shape = alpha, rate = betas[j, ])
    }

    return(apply(mask.dens, 2, prod))
}

## Calculating density due to signal strength
ss.density <- function(fit, session, mask, dists, call_id, animal_id=NULL){
    if (is_animal_model(fit) & is.null(animal_id)) {
      stop("'animal_id' argument must be provided for animal id models")
    }
    
    capt <- get_capt_by_id(fit, call_id, animal_id, session)
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)

    # Grab our estimated parameters, along with other model fitting parameters
    det.pars <- as.list(coef(fit, type="fitted"))
    det.pars$cutoff <- get_ss.opts(fit)$cutoff
    detfn <- get_detfn(fit)
    ss.link <- get_ss_link(fit)
    n.traps <- nrow(get_trap(fit)[[session]])
    
    ss.capt <- capt$ss
    triggered_traps <- capt$trap
    untriggered_traps <- setdiff(seq(n.traps), triggered_traps)
    
    # Initialize density matrix for this call
    mask.dens <- matrix(0, nrow = n.traps, ncol = nrow(mask))
    
    # Similar to how we calculate contribution from detection history alone
    for (i in seq_along(triggered_traps)) {
      mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[triggered_traps[i], ]
      mask.dens[triggered_traps[i], ] <- dnorm(ss.capt[i], mu.ss, det.pars[["sigma.ss"]])
    }
    for (i in seq_along(untriggered_traps)) {
      mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[untriggered_traps[i], ]
      mask.dens[untriggered_traps[i], ] <- pnorm(det.pars$cutoff, mu.ss, det.pars[["sigma.ss"]])
    }
  
    return(apply(mask.dens, 2, prod))
}

## Calculating density due to toa
toa.density <- function(fit, session, mask, dists, call_id, animal_id=NULL){
    if (is_animal_model(fit) & is.null(animal_id)) {
      stop("'animal_id' argument must be provided for animal id models")
    }
  
    capt <- get_capt_by_id(fit, call_id, animal_id, session)
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)

    sigma.toa <- coef(fit, type="fitted")[["sigma.toa"]]

    toa.capt <- capt$toa
    dists.mask <- dists[bincapt == 1, ]
    prod.times <- toa.capt - dists.mask/fit$args$sound.speed
    toa.ssq <- apply(prod.times, 2, function(x) sum((x - mean(x))^2))
    out <- (2*pi*sigma.toa^2)^((1 - sum(bincapt))/2)*
      exp(toa.ssq/(-2*sigma.toa^2))
    
    return(out)
}

## Plots arrows on traps where a detection was made, showing estimated bearing.
estimated_bearing_arrow_plot <- function(fit, session, call_id, animal_id=NULL) {
    if (is_animal_model(fit) & is.null(animal_id)) {
      stop("'animal_id' argument must be provided for animal id models")
    }
    
    capt <- get_capt_by_id(fit, call_id, animal_id, session)
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
    
    is.dist <- !is.null(capt$dist)

    trappos <- get_trap(fit)[[session]]
    
    bearing <- capt$bearing
    dist <- capt$dist
    
    if (is.dist) {
      arrow_len = dist
    } else {
      arrow_len = 0.382 * get_buffer(fit)
    }
    
    activated_traps <- trappos[which(bincapt == 1), , drop = FALSE]
    sinb <- sin(bearing)*arrow_len
    cosb <- cos(bearing)*arrow_len
    
    arrows.df <- data.frame(
      x = activated_traps[, 1],
      y = activated_traps[, 2],
      xend = activated_traps[, 1] + sinb,
      yend = activated_traps[, 2] + cosb
    )
    
    return(geom_segment(data = arrows.df, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                        col="red",
                        arrow = arrow(
                          length = unit(0.2, "cm") # Arrow head
                        )))
}

## Plots circles around traps where a detection was made, showing estimated distance.
estimated_distance_plot <- function(fit, session, call_id, animal_id=NULL){
  if (is_animal_model(fit) & is.null(animal_id)) {
    stop("'animal_id' argument must be provided for animal id models")
  }
  
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)
  
  dist <- capt$dist
  trappos <- get_trap(fit)[[session]]
  
  activated_traps <- trappos[which(bincapt == 1), , drop = FALSE]
  circles.df <- data.frame(
    x = activated_traps[, 1], y = activated_traps[, 2], size=dist
  )
  
  return(geom_point(
    data = circles.df,
    aes(x = x, y = y),
    shape = 21,
    color = "red",
    size = circles.df$size,  
    stroke = 1 
  ))
}

circles <- function(centres, radius, ...){
    bearings <- seq(0, 2*pi, length.out = 100)
    xs <- centres[1] + sin(bearings)*radius
    ys <- centres[2] + cos(bearings)*radius
    lines(xs, ys, ...)
}


