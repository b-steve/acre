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
                      keep.densities = FALSE) {
  
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
    capt.all <- get_capt_for_plot(fit$args)
    capt.all <- capt.all[capt.all$session == session, ]
    
    animal.model <- is_animal_model(fit)
    
    traps <- get_trap(fit)[[session]]
    
    if(animal.model && call_id != "all" && is.null(animal_id)) {
      stop("'animal_id' must be provided when specifying the 'call_id' argument for animal ID models")
    }
    
    if(!is.null(animal_id) && !animal.model) {
      stop("'animal_id' argument should only be provided for animal models. Use 'call_id' argument instead.")
    }
    
    if (!is.null(animal_id)) {
      capt.all <- subset(capt.all, capt.all$animal_ID == animal_id)
    }
    
    ## Setting id properly if "all" selected.
    if (any(call_id == "all")){
      if (animal.model) {
        if(!is.null(animal_id)) {
          call_id <- unique(capt.all$ID)
        } else {
          call_id <- unique(capt.all$animal_ID)
        }
      } else {
        call_id <- unique(capt.all$ID)
      }
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
    base.plot <- base_plot(capt.all, as.data.frame(mask), xlim, ylim)
    
    
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
    if (missing(cols)){
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
    
    densities.list <- list()
    for (i in call_id){
        
        # Sort out plot title
        if (animal.model) {
          if (is.null(animal_id)) {
            sub_title = paste0("session: ", session, ", animal ID: ", i)
          } else {
            sub_title = paste0("session: ", session, ", animal ID: ", animal_id,
                               ", call ID: ", i)
          }
        } else {
          sub_title = paste0("session: ", session, ", call ID: ", i)
        }
        
        base.plot <- base.plot + labs(subtitle = sub_title)
      
        densities <- data.frame(x = mask[, 1], y = mask[, 2])
        
        if (plot.types["combined"]){
            if ((!combine) | (combine & plot.types["capt"])){
                f.combined <- f.x
            } else {
                f.combined <- 0*f.x + 1
            }
        }
        
        # If it's an animal model, but we want to loop over calls separately
        if (!is.null(animal_id)) {
          bincapt <- get_bincapt_by_id(fit, animal_id, session)
          bincapt <- matrix(bincapt[i, ], nrow=1)
        } else {
          # Otherwise just grab it normally
          bincapt <- get_bincapt_by_id(fit, i, session)
        }
        
        # vector of length n.traps
        # 1 if trap was triggered for this id, 0 if not
        traps_triggered <- apply(bincapt, 2, function(col) as.integer(any(col == 1)))
        
        ## Contour due to capture history.
        if (plot.types["capt"] | plot.types["combined"] | plot.types["ss"]){
            det.probs <- det_prob(detfn, det.pars, dists, ss.link)
            
            if(detfn == 'ss'){
              det.probs = 1 - pnorm(det.pars$cutoff, mean = det.probs, sd = det.pars$sigma.ss)
            }
            
            f.capt <- 1
            
            # Note: if animal.model, then bincapt will have multiple rows
            for (row in 1:nrow(bincapt)) {
              # for each capture history specific to this id
              # calculate the contribution to the capture density
              f.capt_unit_contribution <- apply((det.probs ^ bincapt[row, ]) * ((1 - det.probs) ^ (1 - bincapt[row, ])), 2, prod)
              f.capt <- f.capt * f.capt_unit_contribution
            }
            densities$capt <- f.capt
            
            # Previously was 
            # f.capt <- colProds(det.probs*capt + (1 - det.probs)*(1 - capt))
            # f.capt <- apply((det.probs ^ capt) * ((1 - det.probs) ^ (1 - capt)), 2, prod)
      
            if (plot.types["capt"]){
                if (!combine){
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
                # Cool tech here, ifelse() can't return NULL, switch() saves the day
                # is.null(animal_id) + 1 evaluates to 1 if null, 2 if not.
                f.ss.capt <- ss.density(fit = fit, 
                                        id = switch(is.null(animal_id) + 1, animal_id, i), 
                                        session = session, 
                                        mask = mask, 
                                        dists = dists, 
                                        call_id = switch(is.null(animal_id) + 1, i, NULL))
                
                f.ss <- f.ss.capt/f.capt
                densities$ss <- t(f.ss)
                ## Such a hack, but this keeps f.combined correct,
                ## below.
                f.capt <- f.ss.capt
                ## Ughhhh sorry about this one.
                f.ss[f.ss == Inf] <- 0
                if (plot.types["ss"]){
                    if (!combine){
                      contours <- calculate.contour(mask = mask, dens = f.x*f.ss, levels = levels,
                                     nlevels = nlevels, prob = !density)
                      
                      ss_normalised_density <- t(f.x*f.ss / (attr(mask, "a")*10000*sum(f.x*f.ss)))
                      
                      ss_contour <- ggplot2::geom_contour(
                        data = data.frame(x=mask[,1], y=mask[,2], z=ss_normalised_density),
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
            } else if (plot.types["combined"] & combine){
                f.true.capt <- f.capt/f.ss
                f.true.capt[f.ss == 0] <- 0
                if (plot.types["capt"] & plot.types["ss"]){
                    f.combined <- f.combined*f.capt
                } else if (plot.types["capt"] & !plot.types["ss"]){
                    f.combined <- f.combined*f.true.capt
                } else if (!plot.types["capt"] & plot.types["ss"]){
                    f.combined <- f.combined*f.ss
                }
            }
        }
        ## Contour due to estimated bearings.
        if (plot.types["bearing"] | plot.types["combined"] & fit$fit.types["bearing"]){
            f.bearing <- bearing.density(fit = fit, 
                                         id = switch(is.null(animal_id) + 1, animal_id, i), 
                                         session = session, 
                                         mask = mask, 
                                         call_id = switch(is.null(animal_id) + 1, i, NULL))
            densities$bearing <- t(f.bearing)
            if (plot.types["bearing"]){
                if (!combine){
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
                }
            }
        }
        ## Contour due to estimated distances.
        if (plot.types["dist"] | plot.types["combined"] & fit$fit.types["dist"]){
            f.dist <- dist.density(fit = fit, 
                                   id = switch(is.null(animal_id) + 1, animal_id, i), 
                                   session = session, 
                                   mask = mask, 
                                   dists = dists, 
                                   call_id = switch(is.null(animal_id) + 1, i, NULL))
            
            densities$dist <- t(f.dist)
            if (plot.types["dist"]){
                if (!combine){
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
                }
            }
        }
        ## Contour due to measured times of arrival.
        # Note toa data not useful if only 1 trap is triggered
        if ((plot.types["toa"] | plot.types["combined"]) &
            fit$fit.types["toa"] & sum(traps_triggered) > 1){
            f.toa <- toa.density(fit, i, session, mask, dists)
            densities$toa <- f.toa
            if (plot.types["toa"]){
                if (!combine){
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
                  }
              }
        }
        ## Combined contour.
        if (plot.types["combined"]){
            densities$combined <- t(f.combined)
            
            contours <- calculate.contour(mask = mask, dens = f.combined, levels = levels,
                         nlevels = nlevels, prob = !density)
            
            combined_normalised_density <- t(f.combined / (attr(mask, "a")*10000*sum(f.combined))) 
            
            combined_contour <- ggplot2::geom_contour(
              data = data.frame(x=mask[,1], y=mask[,2], z=combined_normalised_density),
              aes(x = x, y = y, z = z),
              breaks=contours$levels, 
              color = cols$combined,
              linetype= ltys$combined
            )
            
            base.plot <- base.plot + combined_contour
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
        
        # Keep track of densities
        densities.list[[as.character(i)]] <- densities
        
        ## Plotting traps, and circles around them.
        if (!add){
          if (plot.ss) {
            plot_traps <- data.frame(traps[which(bincapt == 0), , drop = FALSE])
          } else {
            plot_traps <- data.frame(traps)
          }
          # Plot traps
          base.plot = base.plot + 
            geom_point(data = plot_traps, mapping = aes(x = x, y = y, shape = 4), 
                       size = 3, col="red", stroke=1) + 
            scale_shape_identity()
          
          if (circle.traps & !plot.ss){
            base.plot = base.plot + geom_point(
              data = subset(traps, traps_triggered == 1),
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
                                                                   id = switch(is.null(animal_id) + 1, animal_id, i),
                                                                   session = session,
                                                                   call_id = switch(is.null(animal_id) + 1, i, NULL))
            }
          }
          
          
          # Plotting arrows for estimated bearings.
          if (fit$fit.types["dist"]) {
            if (plot.circles & !plot.arrows){
              base.plot = base.plot + estimated_distance_plot(fit = fit, 
                                                              id = switch(is.null(animal_id) + 1, animal_id, i),
                                                              session = session,
                                                              call_id = switch(is.null(animal_id) + 1, i, NULL))
            }
          }
          
          # Plotting arrows for estimated bearings.
          if (fit$fit.types["ss"]) {
            if (plot.ss){
              # activated_traps <- traps[which(bincapt == 1), , drop = FALSE]
              # base.plot = base.plot + estimated_ss_plot(fit = fit,
              #                                           id = switch(is.null(animal_id) + 1, animal_id, i),
              #                                           session = session,
              #                                           call_id = switch(is.null(animal_id) + 1, i, NULL))
              ss_circles <- estimated_ss_plot(fit = fit,
                                              id = switch(is.null(animal_id) + 1, animal_id, i),
                                              session = session,
                                              call_id = switch(is.null(animal_id) + 1, i, NULL))

              for (circle in ss_circles) {
                base.plot = base.plot + circle
              }

              base.plot = base.plot +
                ggplot2::scale_fill_identity() +                       # Use the fill_color as is
                ggplot2::scale_alpha(range = c(0.1, 0.6)) + # Control alpha range and remove legend
                ggplot2::theme(
                  legend.position = "none"                   # Remove legend if not needed
                )
            }
          }
        }
        
        
        
        plot(base.plot)
        
        if(!add) {
          # Reset the plot
          base.plot = base_plot(capt.all, as.data.frame(mask), xlim, ylim)
        }
        
        if (ask) {
          # Make sure to only ask if we are plotting more than one plot, 
          # and it is not the last plot in the list
          if (length(call_id) > 1 && i != call_id[[length(call_id)]]) {
            prompt_user_for_next_plot()
          }
        }
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
    if (keep.densities) {
      out$densities <- densities.list
    }
    
    else {
      out <- TRUE
    }
    
    invisible(out)
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
bearing.density <- function(fit, id, session, mask, call_id=NULL){
    capt <- get_capt_by_id(fit, id, session)
    bincapt <- get_bincapt_by_id(fit, id, session)
    
    if (!is.null(call_id)) {
      if (is_animal_model(fit)) {
        capt <- subset(capt, capt$ID == call_id)
        bincapt <- matrix(bincapt[call_id, ], nrow=1)
      }
      else {
        stop("'call_id' argument should only be provided for animal models. Use the 'id' argument instead.")
      }
    }
    
    # Initialize final density matrix
    dens <- matrix(1, nrow = 1, ncol = nrow(mask))
    
    # For non-animal_ID models, there is only 1 call
    for (i in 1:length(unique(capt$ID))) {
      call_ID <- unique(capt$ID)[i]
      
      bearing.capt <- subset(capt, capt$ID == call_ID)$bearing
      
      # Initialize density matrix for this call
      mask.dens <- matrix(0, nrow = sum(bincapt[i, ]), ncol = nrow(mask))
      
      # Grab our estimated kappa parameter
      kappa <- coef(fit, type="fitted")[["kappa"]]
      
      # Calculate bearings for each triggered trap
      mask.bearings <- bearings(get_trap(fit)[[session]][bincapt[i, ] == 1, , drop = FALSE], as.matrix(mask))
      
      # Calculate the density for each of our detections
      for (j in 1:sum(bincapt[i, ])) {
        mask.dens[j, ] <- CircStats::dvm(bearing.capt[j], mu = mask.bearings[j, ], kappa = kappa)
      }
      
      # Update the final density matrix
      dens <- dens * apply(mask.dens, 2, prod)
    }
    
    return(dens)
}

## Calculating density due to estimated distances.
dist.density <- function(fit, id, session, mask, dists, call_id=NULL){
    capt <- get_capt_by_id(fit, id, session)
    bincapt <- get_bincapt_by_id(fit, id, session)
    
    if (!is.null(call_id)) {
      if (is_animal_model(fit)) {
        capt <- subset(capt, capt$ID == call_id)
        bincapt <- matrix(bincapt[call_id, ], nrow=1)
      }
      else {
        stop("'call_id' argument should only be provided for animal models. Use the 'id' argument instead.")
      }
    }

    # Initialize final density matrix
    dens <- matrix(1, nrow = 1, ncol = nrow(mask))

    # For non-animal_ID models, there is only 1 call
    for (i in 1:length(unique(capt$ID))) {
      call_ID <- unique(capt$ID)[i]
      dist.capt <- subset(capt, capt$ID == call_ID)$dist

      mask.dens <- matrix(0, nrow = sum(bincapt[i, ]), ncol = nrow(mask))
      
      # Grab our estimated alpha parameter
      alpha <- coef(fit, type="fitted")[["alpha"]]
      betas <- alpha/dists
      
      # Calculate bearings for each triggered trap
      for (j in 1:sum(bincapt[i, ])) {
        mask.dens[j, ] <- dgamma(dist.capt[j], shape = alpha, rate = betas[j, ])
      }

      # Update the final density matrix
      dens <- dens * apply(mask.dens, 2, prod)
    }

    return(dens)
}

## Calculating density due to signal strength
ss.density <- function(fit, id, session, mask, dists, call_id=NULL){
    capt <- get_capt_by_id(fit, id, session)
    bincapt <- get_bincapt_by_id(fit, id, session)
    
    if (!is.null(call_id)) {
      if (is_animal_model(fit)) {
        capt <- subset(capt, capt$ID == call_id)
        # TODO: Error here when call_id specified is no good 
        bincapt <- matrix(bincapt[call_id, ], nrow=1)
      }
      else {
        stop("'call_id' argument should only be provided for animal models. Use the 'id' argument instead.")
      }
    }

    # Initialize final density matrix
    dens <- matrix(1, nrow = 1, ncol = nrow(mask))

    # Grab our estimated parameters, along with other model fitting parameters
    det.pars <- as.list(coef(fit, type="fitted"))
    det.pars$cutoff <- get_ss.opts(fit)$cutoff
    detfn <- get_detfn(fit)
    ss.link <- get_ss_link(fit)
    n.traps <- nrow(get_trap(fit)[[session]])

    for (call_ID in unique(capt$ID)) {
      ss.capt <- subset(capt, capt$ID == call_ID)$ss
      triggered_traps <- subset(capt, capt$ID == call_ID)$trap
      untriggered_traps <- setdiff(seq(n.traps), triggered_traps)

      mask.dens <- matrix(0, nrow = n.traps, ncol = nrow(mask))
      
      # Similar to how we calculate contribution from detection history alone
      for (i in seq_along(triggered_traps)) {
        mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[triggered_traps[i], ]
        mask.dens[triggered_traps[i], ] <- dnorm(ss.capt[i], mu.ss, det.pars[["sigma.ss"]])
      }
      for (i in seq_along(untriggered_traps)) {
        mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[untriggered_traps[i], ]
        mask.dens[untriggered_traps[i], ] <- pnorm(det.pars$cutoff, mu.ss, det.pars[["sigma.ss"]])
        # mask.dens[untriggered_traps[i], ] <- 1 - det_prob(detfn, det.pars, dists[untriggered_traps[i], ], ss.link)
      }

      # Update the final density matrix
      dens <- dens * apply(mask.dens, 2, prod)
    }

    return(dens)
}

## Calculating density due to toa
toa.density <- function(fit, id, session, mask, dists){
    capt <- get_capt_by_id(fit, id, session)
    bincapt <- get_bincapt_by_id(fit, id, session)

    sigma.toa <- coef(fit, type="fitted")[["sigma.toa"]]

    out <- 1
    for (i in 1:length(unique(capt$ID))) {
      call_ID <- unique(capt$ID)[i]
      toa.capt <- subset(capt, capt$ID == call_ID)$toa
      mask.dens <- matrix(0, nrow = sum(bincapt[i, ]), ncol = nrow(mask))

      dists.mask <- dists[bincapt[i, ] == 1, ]

      prod.times <- toa.capt - dists.mask/fit$args$sound.speed
      toa.ssq <- apply(prod.times, 2, function(x) sum((x - mean(x))^2))
      out <- out * (2*pi*sigma.toa^2)^((1 - sum(bincapt))/2)*
        exp(toa.ssq/(-2*sigma.toa^2))
    }
    
    return(out)
}

## Plots arrows on traps where a detection was made, showing estimated bearing.
estimated_bearing_arrow_plot <- function(fit, id, session, call_id=NULL){
    capt <- get_capt_by_id(fit, id, session)
    bincapt <- get_bincapt_by_id(fit, id, session)
    
    if (!is.null(call_id)) {
      if (is_animal_model(fit)) {
        capt <- subset(capt, capt$ID == call_id)
        # TODO: Error here when call_id specified is no good 
        bincapt <- matrix(bincapt[call_id, ], nrow=1)
      }
      else {
        stop("'call_id' argument should only be provided for animal models. Use the 'id' argument instead.")
      }
    }
    
    is.dist <- !is.null(capt$dist)

    trappos <- get_trap(fit)[[session]]
    
    arrows.df <- data.frame()
    for (i in 1:length(unique(capt$ID))) {
      call_ID <- unique(capt$ID)[i]
      bearing <- subset(capt, capt$ID == call_ID)$bearing
      dist <- subset(capt, capt$ID == call_ID)$dist
      
      if (is.dist) {
        arrow_len = dist
      } else {
        arrow_len = 0.382 * get_buffer(fit)
      }

      activated_traps <- trappos[which(bincapt[i, ] == 1), , drop = FALSE]
      sinb <- sin(bearing)*arrow_len
      cosb <- cos(bearing)*arrow_len
      
      arrows.df <- rbind(arrows.df,
                         data.frame(
                           x = activated_traps[, 1],
                           y = activated_traps[, 2],
                           xend = activated_traps[, 1] + sinb,
                           yend = activated_traps[, 2] + cosb
                         ))
    }
    
    return(geom_segment(data = arrows.df, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                        col="red",
                        arrow = arrow(
                          length = unit(0.2, "cm") # Arrow head
                        )))
}

## Plots circles around traps where a detection was made, showing estimated distance.
estimated_distance_plot <- function(fit, id, session, call_id=NULL){
  capt <- get_capt_by_id(fit, id, session)
  bincapt <- get_bincapt_by_id(fit, id, session)
  
  if (!is.null(call_id)) {
    if (is_animal_model(fit)) {
      capt <- subset(capt, capt$ID == call_id)
      # TODO: Error here when call_id specified is no good 
      bincapt <- matrix(bincapt[call_id, ], nrow=1)
    }
    else {
      stop("'call_id' argument should only be provided for animal models. Use the 'id' argument instead.")
    }
  }
  
  circles.df <- data.frame()
  for (i in 1:length(unique(capt$ID))) {
    call_ID <- unique(capt$ID)[i]
    dist <- subset(capt, capt$ID == call_ID)$dist

    trappos <- get_trap(fit)[[session]]
    
    activated_traps <- trappos[which(bincapt[i, ] == 1), , drop = FALSE]
    
    circles.df <- rbind(circles.df,
                        data.frame(
                          x = activated_traps[, 1], y = activated_traps[, 2], size=dist
                        ))
  }
  
  return(geom_point(
    data = circles.df,
    aes(x = x, y = y),
    shape = 21,
    color = "red",
    size = circles.df$size,  
    stroke = 1 
  ))
}

estimated_ss_plot <- function(fit, id, session, call_id=NULL) {
  capt <- get_capt_by_id(fit, id, session)
  bincapt <- get_bincapt_by_id(fit, id, session)
  
  if (!is.null(call_id)) {
    if (is_animal_model(fit)) {
      capt <- subset(capt, capt$ID == call_id)
      # TODO: Error here when call_id specified is no good 
      bincapt <- matrix(bincapt[call_id, ], nrow=1)
    }
    else {
      stop("'call_id' argument should only be provided for animal models. Use the 'id' argument instead.")
    }
  }
  
  circles.list <- list()
  for (i in 1:length(unique(capt$ID))) {
    call_ID <- unique(capt$ID)[i]
    ss <- subset(capt, capt$ID == call_ID)$ss
    
    trappos <- get_trap(fit)[[session]]
    
    activated_traps <- trappos[which(bincapt[i, ] == 1), , drop = FALSE]
    
    
    # return(geom_point(data = activated_traps, mapping = aes(x = activated_traps[,1], 
    #                                                  y = activated_traps[,2], 
    #                                                  colour = capt$ss), 
    #            size = 5))
    # rbind(circles.df,
    #       data.frame(
    #         x = activated_traps[, 1], y = activated_traps[, 2], size=dist
    #       ))
    for (j in 1:nrow(activated_traps)) {
      trap_x <- activated_traps[j, 1]
      trap_y <- activated_traps[j, 2]

      circle_data <- generate_concentric_circles(
        num_circles = 4,          # More circles for a smoother gradient
        max_radius = 1.5,            # Adjust as needed
        min_radius = 0.5,
        center_x = trap_x,              # Center of the circles
        center_y = trap_y,
        points_per_circle = 50,   # Higher for smoother circles
        fill_color = "red",
        alpha_start = 0.05,
        alpha_end = 0.5
      )

      circles.list[[j]] <- ggplot2::geom_polygon(
        data = circle_data,
        aes(x = x, y = y, group = group, fill = fill_color, alpha = alpha),
        color = NA
      )
    }
  }
  
  return(circles.list)
  
}

circles <- function(centres, radius, ...){
    bearings <- seq(0, 2*pi, length.out = 100)
    xs <- centres[1] + sin(bearings)*radius
    ys <- centres[2] + cos(bearings)*radius
    lines(xs, ys, ...)
}

# Function to generate concentric circles data
generate_concentric_circles <- function(
    num_circles = 4,          # Number of concentric circles
    max_radius = 1,            # Radius of the largest circle
    min_radius = 0.1,         # Radius of the smallest circle
    center_x = 0,              # X-coordinate of the center
    center_y = 0,              # Y-coordinate of the center
    points_per_circle = 50,   # Number of points to define each circle
    fill_color = "red",    # Fill color for the circles
    alpha_start = 0.1,         # Alpha for the largest circle
    alpha_end = 1               # Alpha for the smallest circle
) {
  
  # Function to generate a single circle's coordinates
  generate_circle <- function(cx, cy, radius, points) {
    angles <- seq(0, 2 * pi, length.out = points)
    data.frame(
      x = cx + radius * cos(angles),
      y = cy + radius * sin(angles)
    )
  }
  
  # Sequence of radii from max_radius to min_radius
  radii <- seq(max_radius, min_radius, length.out = num_circles)
  
  # Corresponding alpha values from alpha_start to alpha_end
  alphas <- seq(alpha_start, alpha_end, length.out = num_circles)
  
  # Initialize an empty data frame to store all circles
  all_circles <- data.frame()
  
  # Generate data for each circle
  for (i in seq_along(radii)) {
    circle <- generate_circle(center_x, center_y, radii[i], points_per_circle)
    circle$group <- i               # Group identifier for each circle
    circle$alpha <- alphas[i]       # Alpha value for each circle
    circle$fill_color <- fill_color # Fill color
    all_circles <- rbind(all_circles, circle)
  }
  
  # Return the complete data frame
  return(all_circles)
}
