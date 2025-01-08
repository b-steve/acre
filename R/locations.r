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
locations <- function(fit, id = "all", session = 1, infotypes = NULL,
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
                      plot.circles = "dist" %in% fit$infotypes &
                          !("bearing" %in% fit$infotypes),
                      arrow.length = NULL, show.legend = FALSE,
                      show.axes = TRUE, add = FALSE, ask = TRUE,
                      keep.densities = FALSE){
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
    if (missing(infotypes) & combine){
        stop("Argument `combine' is only useful if `infotypes' is provided.")
    }
    ## Error if new mask covariates are not provided when they need to be.
    is.new.mask <- !missing(mask)
    if (fit$fit.ihd & is.new.mask & is.null(newdata)){
        stop("Covariate values for the mask object must be provided via the `newdata' argument.")
    }
    
    ## Extracting the session's capture history.
    capt.all <- get_capt(fit, session)
    animal.model <- is_animal_model(fit)
    
    ## Setting id properly if "all" selected.
    if (id == "all"){
      if (animal.model) {
        id <- unique(capt.all$animal_ID)
      } else {
        id <- 1:nrow(capt.all$bincapt)
      }
    }

    ## Saving estimated locations.
    if (keep.estlocs){
        estlocs <- matrix(0, nrow = length(id), ncol = 2)
        j <- 1
    }
    
    ## Sorting out limits.
    if (is.null(xlim)){
        xlim <- range(mask[,1])
    }
    if (is.null(ylim)){
        ylim <- range(mask[,2])
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
    traps <- get_trap(fit)[[session]]
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
    for (i in id){
        densities <- data.frame(x = mask[, 1], y = mask[, 2])
        
        if (plot.types["combined"]){
            if ((!combine) | (combine & plot.types["capt"])){
                f.combined <- f.x
            } else {
                f.combined <- 0*f.x + 1
            }
        }
        
        bincapt <- get_bincapt_by_id(fit, i)
        
        # vector of length n.traps
        # 1 if trap was triggered for this id, 0 if not
        traps_triggered <- apply(bincapt, 2, function(col) as.integer(any(col == 1)))
        
        ## Contour due to capture history.
        if (plot.types["capt"] | plot.types["combined"] | plot.types["ss"]){
            det.probs <- det_prob(detfn, det.pars, dists, ss.link)
            
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
                    show.contour(mask = mask, dens = f.x*f.capt, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$capt,
                                 lty = ltys$capt, show.labels = show.labels,
                                 plot.contours = plot.contours)
                }
            }
            if (fit$fit.types["ss"]){
                f.ss.capt <- ss.density(fit, i, session, mask, dists)
                f.ss <- f.ss.capt/f.capt
                densities$ss <- t(f.ss)
                ## Such a hack, but this keeps f.combined correct,
                ## below.
                f.capt <- f.ss.capt
                ## Ughhhh sorry about this one.
                f.ss[f.ss == Inf] <- 0
                if (plot.types["ss"]){
                    if (!combine){
                        show.contour(mask = mask, dens = f.x*f.ss, levels = levels,
                                     nlevels = nlevels, prob = !density, col = cols$ss,
                                     lty = ltys$ss, show.labels = show.labels,
                                     plot.contours = plot.contours)
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
            f.bearing <- bearing.density(fit, i, session, mask)
            densities$bearing <- t(f.bearing)
            if (plot.types["bearing"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.bearing, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$bearing,
                                 lty = ltys$bearing, show.labels = show.labels,
                                 plot.contours = plot.contours)
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
            f.dist <- dist.density(fit, i, session, mask, dists)
            densities$dist <- t(f.dist)
            if (plot.types["dist"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.dist, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$dist,
                                 lty = ltys$dist, show.labels = show.labels,
                                 plot.contours = plot.contours)
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
        if (plot.types["toa"] | plot.types["combined"] &
            fit$fit.types["toa"] & sum(traps_triggered) > 1){
            f.toa <- toa.density(fit, i, session, mask, dists)
            densities$toa <- f.toa
            if (plot.types["toa"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.toa, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$toa,
                                 lty = ltys$toa, show.labels = show.labels,
                                 plot.contours = plot.contours)
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
            show.contour(mask = mask, dens = f.combined, levels = levels,
                         nlevels = nlevels, prob = !density, col = cols$combined,
                         lty = ltys$combined, show.labels = show.labels,
                         plot.contours = plot.contours)
        }
        if (plot.types["mrds"]){
            loc <- capt.all$mrds[id, , drop = FALSE]
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
          points(traps, col = trap.col, pch = 4, lwd = 2)
          if (circle.traps){
            # TODO: Figure out why this condition was needed
            # if (length(id) == 1){
            # points(traps[traps_triggered == 1, , drop = FALSE], col = trap.col, cex = 2, lwd = 2)
            # }
            points(traps[traps_triggered == 1, , drop = FALSE], col = trap.col, cex = 2, lwd = 2)
          }
        }
        ## Plotting arrows for estimated bearings.
        if (fit$fit.types["bearing"]){
          if (plot.arrows){
            if (fit$fit.types["dist"]){
              if (animal.model) {
                arrow.length <- as.matrix(
                  reshape(
                  subset(capt.all, animal_ID == i)[, c("ID", "trap", "dist")], timevar = "trap", 
                                                  idvar = "ID", direction = "wide")[, -1])
                
                } else {
                arrow.length <- matrix(capt.all$dist[i, ], nrow=1)
              }
            }
            show.arrows(fit, i, session, arrow.length, trap.col)
          }
        }
        # Plotting circles for estimated distances.
        if (fit$fit.types["dist"]){
          if (plot.circles & !plot.arrows){
            show.circles(fit, i, session, trap.col)
          }
        }
        
        if (ask) {
          # Make sure to only ask if we are plotting more than one plot, 
          # and it is not the last plot in the list
          if (length(id) > 1 && i != id[[length(id)]]) {
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

## Helper to get stuff in the right form for contour().
#' Title
#'
#' @export
show.contour <- function(mask, dens, nlevels, levels, prob, col = "black", lty = 1, show.labels, plot.contours,
                         return.contours = FALSE){
    if (plot.contours){
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
        contour(x = unique.x, y = unique.y, z = z, levels = levels, labels = labels,
                col = col, lty = lty, drawlabels = show.labels, add = F)
    }
}

## Calculating density due to estimated bearings.
bearing.density <- function(fit, id, session, mask){
    animal.model <- is_animal_model(fit)
    
    capt <- get_capt(fit, session)
    n.calls <- get_n_calls(fit, id, session)
    
    if (animal.model) {
      capt <- subset(capt, capt$animal_ID == id)
    }
    
    dens <- matrix(1, nrow = 1, ncol = nrow(mask))
    
    for (call in 1:n.calls) {
      if (animal.model) {
        id <- unique(capt$ID)[[call]]
        bincapt <- subset(capt, capt$ID == id)$bincapt
        bearing.capt <- subset(capt, capt$ID == id)$bearing[bincapt == 1]
      } else {
        bincapt <- capt$bincapt[id, ]
        bearing.capt <- capt$bearing[id, bincapt == 1]
      }
      
      mask.dens <- matrix(0, nrow = sum(bincapt), ncol = nrow(mask))
      
      kappa <- coef(fit, type="fitted")[["kappa"]]
      mask.bearings <- bearings(get_trap(fit)[[session]][bincapt == 1, , drop = FALSE], mask)
      for (i in 1:sum(bincapt)) {
        mask.dens[i, ] <- CircStats::dvm(bearing.capt[i], mu = mask.bearings[i, ], kappa = kappa)
      }
      
      dens <- dens * apply(mask.dens, 2, prod)
    }
    
    return(dens)
}

## Calculating density due to estimated distances.
dist.density <- function(fit, id, session, mask, dists){
    animal.model <- is_animal_model(fit)

    capt <- get_capt(fit, session)
    n.calls <- get_n_calls(fit, id, session)

    if (animal.model) {
      capt <- subset(capt, capt$animal_ID == id)
    }

    dens <- matrix(1, nrow = 1, ncol = nrow(mask))

    for (call in 1:n.calls) {
      if (animal.model) {
        id <- unique(capt$ID)[[call]]
        bincapt <- subset(capt, capt$ID == id)$bincapt
        dist.capt <- subset(capt, capt$ID == id)$dist[bincapt == 1]
      } else {
        bincapt <- capt$bincapt[id, ]
        dist.capt <- capt$dist[id, bincapt == 1]
      }

      mask.dens <- matrix(0, nrow = sum(bincapt), ncol = nrow(mask))

      alpha <- coef(fit, type="fitted")[["alpha"]]
      betas <- alpha/dists
      for (i in 1:sum(bincapt)) {
        mask.dens[i, ] <- dgamma(dist.capt[i], shape = alpha, rate = betas[i, ])
      }

      dens <- dens * apply(mask.dens, 2, prod)
    }

    return(dens)
}

ss.density <- function(fit, id, session, mask, dists){
    animal.model <- is_animal_model(fit)

    capt <- get_capt(fit, session)
    n.calls <- get_n_calls(fit, id, session)

    if (animal.model) {
      capt <- subset(capt, capt$animal_ID == id)
    }

    dens <- matrix(1, nrow = 1, ncol = nrow(mask))

    det.pars <- as.list(coef(fit, type="fitted"))
    det.pars$cutoff <- get_ss.opts(fit)$cutoff
    detfn <- get_detfn(fit)
    ss.link <- get_ss_link(fit)
    n.traps <- nrow(get_trap(fit)[[session]])

    for (call in 1:n.calls) {
      if (animal.model) {
        id <- unique(capt$ID)[[call]]
        bincapt <- subset(capt, capt$ID == id)$bincapt
        ss.capt <- subset(capt, capt$ID == id)$ss
      } else {
        bincapt <- capt$bincapt[id, ]
        ss.capt <- capt$ss[id, ]
      }

      mask.dens <- matrix(0, nrow = n.traps, ncol = nrow(mask))
      for (i in 1:n.traps) {
        if (bincapt[i] == 0) {
          # mask.dens[i, ] <- 1 - calc.detfn(dists[i, ], detfn, det.pars, ss.link)
          mask.dens[i, ] <- 1 - det_prob(detfn, det.pars, dists[i, ], ss.link)
        } else if (bincapt[i] == 1) {
          mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[i, ]
          mask.dens[i, ] <- dnorm(ss.capt[i], mu.ss, det.pars[["sigma.ss"]])
        } else {
          stop("The binary capture history must only contain 0s and 1s.")
        }
      }

      dens <- dens * apply(mask.dens, 2, prod)
    }

    return(dens)
}

toa.density <- function(fit, id, session, mask, dists){
    animal.model <- is_animal_model(fit)

    capt <- get_capt(fit, session)
    n.calls <- get_n_calls(fit, id, session)

    if (animal.model) {
      capt <- subset(capt, capt$animal_ID == id)
    }

    sigma.toa <- coef(fit, type="fitted")[["sigma.toa"]]

    out <- 1
    for (call in 1:n.calls) {
      if (animal.model) {
        id <- unique(capt$ID)[[call]]
        bincapt <- subset(capt, capt$ID == id)$bincapt
        toa.capt <- subset(capt, capt$ID == id)$toa[bincapt == 1]
        
        # If only 1 detection, toa doesn't give us any info, so skip
        if (sum(bincapt) <= 1) {
          next 
        }
      } else {
        bincapt <- capt$bincapt[id, ]
        toa.capt <- capt$toa[id, bincapt == 1]
      }
    
      dists.mask <- dists[bincapt == 1, ]

      prod.times <- toa.capt - dists.mask/fit$args$sound.speed
      toa.ssq <- apply(prod.times, 2, function(x) sum((x - mean(x))^2))
      out <- out * (2*pi*sigma.toa^2)^((1 - sum(bincapt))/2)*
        exp(toa.ssq/(-2*sigma.toa^2))
    }
    
    return(out)
}

## Plots arrows on traps where a detection was made, showing estimated bearing.
show.arrows <- function(fit = NULL, id, session, arrow.length = NULL, trap.col, capt = NULL, traps = NULL){
    xlim <- par("usr")[c(1, 2)]
    ylim <- par("usr")[c(3, 4)]

    if (!is.null(fit)){
        animal.model <- is_animal_model(fit)
        capt.all <- get_capt(fit, session)
        bincapt <- get_bincapt_by_id(fit, id, session)
    } else {
        capt.all <- capt
        bincapt <- matrix(capt.all$bincapt[id, ], nrow = 1)
        animal.model <- "animal_ID" %in% colnames(capt.all)
    }
    
    if (animal.model) {
      capt.all <- subset(capt.all, animal_ID == id)
    }
    
    if (animal.model) {
      # Reshape the capt data frame, so each row is a unique CALL, and each 
      # column represents the bearing of the corresponding trap
      # Note that id no bearing was detected, the value will be 0
      # Also we drop the ID column as it is no longer needed
      bearing.capt <- as.matrix(reshape(capt.all[, c("ID", "trap", "bearing")], timevar = "trap", 
              idvar = "ID", direction = "wide")[, -1])
    } else {
      bearing.capt <- matrix(capt.all$bearing[id, ], nrow=1)
    }
    
    for (i in 1:nrow(bearing.capt)) {
      bearing <- bearing.capt[i, bincapt[i, ] == 1]
      
      # TODO: This is definitley not finished
      if (is.null(arrow.length)){
        lengths <- 0.05*min(c(diff(range(xlim)), diff(range(ylim))))
      } else {
        lengths <- arrow.length[i, bincapt[i, ] == 1]
      }
      
      
      if (!is.null(fit)){
        trappos <- get_trap(fit)[[session]]
      } else {
        trappos <- traps
      }
      trappos <- trappos[which(bincapt[i, ] == 1), , drop = FALSE]
      sinb <- sin(bearing)*lengths
      cosb <- cos(bearing)*lengths
      arrows(trappos[, 1], trappos[, 2], trappos[, 1] + sinb, trappos[, 2] + cosb,
             length = 0.1, col = trap.col, lwd = 2)
    }
}

## Plots circles around traps where a detection was made, showing estimated distance.
show.circles <- function(fit = NULL, id, session, trap.col, capt = NULL, traps = NULL){
    if (!is.null(fit)){
      animal.model <- is_animal_model(fit)
      capt.all <- get_capt(fit, session)
      bincapt <- get_bincapt_by_id(fit, id, session)
    } else {
      capt.all <- capt
      bincapt <- matrix(capt.all$bincapt[id, ], nrow = 1)
      animal.model <- "animal_ID" %in% colnames(capt.all)
    }
    
    if (animal.model) {
      capt.all <- subset(capt.all, animal_ID == id)
    }
    
    if (animal.model) {
      # Reshape the capt data frame, so each row is a unique CALL, and each 
      # column represents the bearing of the corresponding trap
      # Note that id no bearing was detected, the value will be 0
      # Also we drop the ID column as it is no longer needed
      dist.capt <- as.matrix(reshape(capt.all[, c("ID", "trap", "dist")], timevar = "trap", 
                                        idvar = "ID", direction = "wide")[, -1])
    } else {
      dist.capt <- matrix(capt.all$dist[id, ], nrow=1)
    }
  
    if (!is.null(fit)){
        trappos <- get_trap(fit)[[session]]
    } else {
        trappos <- traps
    }

    for (i in 1:nrow(bincapt)) {
      positions <- trappos[which(bincapt[i, ] == 1), , drop = FALSE]

      radius <- dist.capt[i, which(bincapt[i, ] == 1)]
      
      for (j in 1:nrow(positions)) {
        circles(positions[j, ], radius[j], col = trap.col, lwd = 2)
      }
    }
}

circles <- function(centres, radius, ...){
    bearings <- seq(0, 2*pi, length.out = 100)
    xs <- centres[1] + sin(bearings)*radius
    ys <- centres[2] + cos(bearings)*radius
    lines(xs, ys, ...)
}
