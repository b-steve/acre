get.capt <- function(fit, session = NULL){
  if (fit$n.sessions == 1){
    out <- fit$args$capt
  } else {
    if (!is.numeric(session)) {
      stop("Must provide session number for multiple session data.")
    }
    out <- subset(fit$args$capt, session == 1)
  }
  out
}

get.mask <- function(fit, session = NULL, as.list = NULL){
  if (length(session) == 1){
    if (session == "all"){
      session <- 1:fit$n.sessions
    }
  }
  if (fit$n.sessions == 1){
    session <- 1
  }
  if (is.null(session) | !is.list(fit$args$mask)){
    out <- fit$args$mask
  } else {
    out <- fit$args$mask[session]
  }
  if (is.null(as.list)){
    as.list <- length(out) > 1
  }
  if (length(session) > 1){
    if (!as.list){
      out <- do.call("rbind", out)
    }
  } else {
    if (!as.list){
      out <- out[[1]]
    }
  }
  out
}

get.traps <- function(fit, session = NULL, as.list = NULL){
  return(as.matrix(fit$args$traps[[1]]))
}




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
                      mask = get.mask(fit, session), newdata = NULL,
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
                      show.axes = TRUE, add = FALSE, ask = TRUE){
    ## Error if session argument is too large.
    if (session > fit$n.sessions){
        if (fit$n.sessions == 1){
            stop(paste("Argument 'session' too large; there was only ", fit$n.sessions, " session.", sep = ""))
        } else {
            stop(paste("Argument 'session' too large; there were only ", fit$n.sessions, " sessions.", sep = ""))
        }
    }
  
    # Check if animal_ID model
    if('animal_ID' %in% colnames(fit$args$capt)){
      animal.model = TRUE
    } else {
      animal.model = FALSE
    }
  
    # JO CHECK
    ## Error for locations() with a directional model.
    if (!is.null(fit$args$ss.opts$directional)){
        if (fit$args$ss.opts$directional){
            stop("The locations() function has not yet been implemented for directional model fits.")
        }
    }
    # JO CHECK (Probably delete this)
    if (!is.null(fit$args$ss.opts$het.source)){
        if (fit$args$ss.opts$het.source){
            stop("The locations() function has not yet been implemented for heterogeneous source strength model fits.")
        }
    }
    # JO CHECK (Do we have this?)
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
    # JO CHECK What is fit.ihd again?
    if (fit$fit.ihd & is.new.mask & is.null(newdata)){
        stop("Covariate values for the mask object must be provided via the `newdata' argument.")
    }
    ## Extracting the session's capture history.
    capt.all <- get.capt(fit, session)
    
    ## Setting id properly if "all" selected.
    if (id == "all"){
        id <- 1:length(capt.all$bincapt) 
    }
    ## Saving estimated locations.
    if (keep.estlocs){
        estlocs <- matrix(0, nrow = length(id), ncol = 2)
        j <- 1
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
        infotypes <- c(fit$infotypes, "capt", "combined"[any.infotypes])
    }
    ## If "ss" is an infotype, set to "capt". OR NOT. <<- Wait what the hell is going on in this comment.
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
    traps <- get.traps(fit, session)
    detfn <- get_detfn(fit)
    ss.link <- get_ss_link(fit)
    if (fit$fit.types["ss"]){
      fit$coefficients[["cutoff"]] <- fit$args$ss.opts$cutoff
    }
    pars = as.list(fit$coefficients)
    
    # dists[i, j] = distance from trap i to mask point j
    dists <- distances(traps, mask)
    
    ## Calculating density due to animal locations.
    det.probs <- det_prob(get_detfn(fit), fit$coefficients, dists, ss.link)
    p.det <- p.dot.defaultD(points = mask, traps = traps, detfn = get_detfn(fit), pars=as.list(fit$coefficients), A=attr(mask, 'area'), ss.link)
    ## Divide by normalising constant; not conversion to square metres.
    a <- attr(mask, "area")
    if (fit$fit.ihd){
        if (is.new.mask){
            D.mask.sess <- predict(fit, newdata)
        } else {
            D.mask.sess <- fit$D.mask[[session]]
        }
        f.x <- p.det*D.mask.sess/(a*sum(p.det*D.mask.sess))
    } else {
        # f.x <- p.det/(a*sum(p.det))
        # the numerator cancels with the f(ωi|xi; γ) term
        if (animal.model) {
          f.x <- 1 / p.det
        } else {
          f.x <- 1 / p.det
        }
        
    }

    # f.x <- conditional_location_density(capt.all, traps, mask, pars, detfn, ss.link)
    
    
    numerator <- function(capt, animal_id) {
      # split capt up 
      animal_capt <- subset(capt, animal_ID == animal_id)
      
      animal_calls <- split(animal_capt, animal_capt$ID)
      
      f.x <- 1
      # for each call
      for (i in 1:length(animal_calls)) {
        call <- animal_calls[i]
        
        # calculate density from capture history
        
        
        
        # calculate density from covariates
      }
      
    }
    
    # numerator()
    
    if(ask){
      ## Setting par(ask).
      ask.save <- par("ask")
      par(ask = TRUE)
      ## Making sure par is restored on function exit.
      on.exit(par(ask = ask.save))
    }
    
    
    ## Calculating conditional density of capture history, given location.
    for (i in id){
      print(i)
        if (plot.types["combined"]){
            if ((!combine) | (combine & plot.types["capt"])){
                f.combined <- f.x
              # f.combined <- 1
            } else {
                f.combined <- 0*f.x + 1
            }
        }
        
        # Sort out co-variate data
        covariates <- list()
        if (fit$fit.types["bearing"]) {
          covariates$bearing <- capt.all$bearing[i, ]
        }
        if (fit$fit.types["dist"]) {
          covariates$dist <- capt.all$dist[i, ]
        }
        if (fit$fit.types["toa"]) {
          covariates$toa <- capt.all$toa[i, ]
        }
        
        sound.speed <- fit$args$sound.speed
      
        # Calculate density
        calculate_densities <- function(covariates, bincapt, traps, mask, pars, detfn, ss.link, sound.speed) {
          # dists[i, j] = distance from trap i to mask point j
          dists <- distances(traps, mask)
          
          # Area of each mask point
          A <- attr(mask, "area")

          location.density <- function(mask, dists, traps, detfn, ss.link, pars, A) {
            # Probability of detection at each mask point
            det.probs <- det_prob(detfn, pars, dists, ss.link)
            
            # Probability detected by at least 1 detector
            p.det <- p.dot.defaultD(mask, traps, detfn, ss.link, pars, A)
            
            # In-homogeneous density magic,
            #       Will come back to fix this once I understand
            if (fit$fit.ihd){
              # if (is.new.mask){
              #   D.mask.sess <- predict(fit, newdata)
              # } else {
              #   D.mask.sess <- fit$D.mask[[session]]
              # }
              # f.x <- p.det*D.mask.sess/(sum(p.det*D.mask.sess))
            } else {
              # the numerator cancels with the f(ωi|xi; γ) term
              f.x <- 1 / p.det
            }
          }
          
          f.location <- location.density(mask, dists, traps, detfn, ss.link, pars, A)
          
          capt.density <- function(bincapt, pars) {
              # f.capt <- apply(det.probs * capt + (1 - det.probs) * (1 - capt), 2, prod)
              return(apply((det.probs ^ bincapt) * ((1 - det.probs) ^ (1 - bincapt)), 2, prod))
          }
          
          f.capt <- capt.density(bincapt, pars)
          # f.ss <- function(fit, i, session, mask, dists) {
          #   f.ss.capt <- ss.density(fit, i, session, mask, dists)
          #   # Calculate density from ss strength alone?
          #   f.ss <- f.ss.capt / f.capt
          #   # Set to appropriate density since we are using ss.
          #   # f.capt <- f.ss.capt
          #   
          #   if (any(f.ss == Inf)) {
          #     stop("Something went wrong calculating ss capture history density (f.ss == infinity).")
          #   }
          #   
          #   return(f.ss)
          # }
          # 
          f.bearing <- ifelse(fit$fit.types["bearing"], bearing.density(bincapt, covariates, traps, mask, pars), 1)
          
          f.toa <- ifelse(fit$fit.types["toa"], toa.density(bincapt, covariates, traps, mask, pars, dists, sound.speed), 1)
          
          f.dist <- ifelse(fit$fit.types["dist"], dist.density(bincapt, covariates, traps, mask, pars, dists), 1)
          
          densities <- list(
            location = f.location,
            capt = f.capt,
            # ss = f.ss,
            bearing = f.bearing,
            toa = f.toa,
            dist = f.dist
          )
          
          return(densities)
        }
        
        # densities <- calculate_densities(covariates, capt.all$bincapt[i,], traps, mask, pars, detfn, ss.link, sound.speed)
        
        # JO CHECK notice that we're assuming different data format for bincapt, make sure this is a valid assumption
        capt <- capt.all$bincapt[i,]
        # capt <- capt.all$bincapt[i]
        
        ## Contour due to capture history.
        if (plot.types["capt"] | plot.types["combined"] | plot.types["ss"]){
            # JO CHECK
            # det.pars <- get.par(fit, fit$detpars, as.list = TRUE)
            det.pars <- fit$coefficients
            if (fit$fit.types["ss"]){
                det.pars$cutoff <- fit$args$ss.opts$cutoff
            }
            # JO CHECK only hn again
            
            f.capt <- apply(det.probs * capt + (1 - det.probs) * (1 - capt), 2, prod)
            # f.capt <- apply((det.probs ^ capt) * ((1 - det.probs) ^ (1 - capt)), 2, prod)
            if (plot.types["capt"]) {
                if (!combine) {
                    show.contour(mask = mask, dens = f.x*f.capt, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$capt,
                                 lty = ltys$capt, show.labels = show.labels,
                                 plot.contours = plot.contours, ask=ask)
                }
            }
            if (fit$fit.types["ss"]){
                f.ss.capt <- ss.density(fit, i, session, mask, dists)
                f.ss <- f.ss.capt/f.capt
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
                                     plot.contours = plot.contours, ask = ask)
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
            # f.bearing <- bearing.density(fit, i, session, mask)
          # bincapt, covariates, traps, mask, pars
            covariates <- data.frame(bearing = capt.all$bearing[i ,])
            f.bearing <- bearing.density(capt.all$bincapt[i,], covariates, traps, mask, pars)
            if (plot.types["bearing"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.bearing, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$bearing,
                                 lty = ltys$bearing, show.labels = show.labels,
                                 plot.contours = plot.contours, ask = ask)
                }
            }
            if (plot.types["combined"]){
                if ((!combine) | (combine & plot.types["bearing"])){
                    f.combined <- f.combined* (f.bearing / sum(f.bearing))
                }
            }
        }
        ## Contour due to estimated distances.
        if (plot.types["dist"] | plot.types["combined"] & fit$fit.types["dist"]){
            covariates <- data.frame(dist = capt.all$dist[i ,])
            # f.dist <- dist.density(fit, i, session, mask, dists)
            f.dist <- dist.density(capt.all$bincapt[i,], covariates, traps, mask, pars, dists)
            if (plot.types["dist"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.dist, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$dist,
                                 lty = ltys$dist, show.labels = show.labels,
                                 plot.contours = plot.contours, ask = ask)
                }
            }
            if (plot.types["combined"]){
                if ((!combine) | (combine & plot.types["dist"])){
                    f.combined <- f.combined*f.dist
                }
            }
        }
        ## Contour due to measured times of arrival.
        if (plot.types["toa"] | plot.types["combined"] &
            fit$fit.types["toa"] & sum(capt) > 1){
            # f.toa <- toa.density(fit, i, session, mask, dists)
            covariates <- data.frame(toa = capt.all$toa[i ,])
            f.toa <- toa.density(capt.all$bincapt[i,], covariates, traps, mask, pars, dists, fit$args$sound.speed)
            if (plot.types["toa"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.toa, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$toa,
                                 lty = ltys$toa, show.labels = show.labels,
                                 plot.contours = plot.contours, ask = ask)
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
            show.contour(mask = mask, dens = f.combined, levels = levels,
                         nlevels = nlevels, prob = !density, col = cols$combined,
                         lty = ltys$combined, show.labels = show.labels,
                         plot.contours = plot.contours, ask = ask)
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
              arrow.length <- capt.all$dist[i, capt == 1]
            }
            show.arrows(fit, i, session, arrow.length, trap.col)
          }
        }
        ## Plotting circles for estimated distances.
        if (fit$fit.types["dist"]){
          if (plot.circles & !plot.arrows){
            show.circles(fit, i, session, trap.col)
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
    if (keep.estlocs){
        out <- estlocs
    } else {
        out <- TRUE
    }
    
    invisible(out)
}

## Helper to get stuff in the right form for contour().
show.contour <- function(mask, dens, nlevels, levels, prob, col = "black", lty = 1, show.labels, plot.contours, ask){
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
                col = col, lty = lty, drawlabels = show.labels, add = !ask)
    }
}

## Calculating density due to estimated bearings.
bearing.density <- function(bincapt, covariates, traps, mask, pars){
  if (!('bearing' %in% names(covariates))) {
    stop("Something went wrong calculating bearing density,
         'bearing' column missing from covariate matrix.")
  }

  bearing.capt <- covariates$bearing[bincapt == 1]
  kappa <- pars$kappa
  mask.bearings <- bearings(traps[bincapt == 1, , drop = FALSE], mask)
  mask.dens <- matrix(0, nrow = sum(bincapt), ncol = nrow(mask))
  for (i in 1:sum(bincapt)){
    mask.dens[i, ] <- CircStats::dvm(bearing.capt[i], mu = mask.bearings[i, ], kappa = kappa)
  }
  ## Returning densities.
  apply(mask.dens, 2, prod)
}

## Calculating density due to estimated distances.
dist.density <- function(bincapt, covariates, traps, mask, pars, dists){
  if (!('dist' %in% names(covariates))) {
    stop("Something went wrong calculating distance density,
         'dist' column missing from covariate matrix.")
  }
  
  dists <- dists[bincapt == 1, , drop = FALSE]
  dist.capt <- covariates$dist[bincapt == 1]
  alpha <- pars$alpha
  mask.dens <- matrix(0, nrow = sum(bincapt), ncol = nrow(mask))
  betas <- alpha / dists
  for (i in 1:sum(bincapt)){
    mask.dens[i, ] <- dgamma(dist.capt[i], shape = alpha, rate = betas[i, ])
  }
  ## Returning densities.
  apply(mask.dens, 2, prod)
}

ss.density <- function(fit, id, session, mask, dists){
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    ss.capt <- capt.all$ss[id, ]
    det.pars <- as.list(fit$coefficients)
    detfn <- get_detfn(fit)
    ss.link <- get_ss_link(fit)
    n.traps <- nrow(get.traps(fit, session))
    mask.dens <- matrix(0, nrow = n.traps, ncol = nrow(mask))
    for (i in 1:n.traps){
        if (capt[i] == 0){
            mask.dens[i, ] <- 1 - det_prob(detfn, det.pars,dists[i, ], ss.link=ss.link)
        } else if (capt[i] == 1){
            mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[i, ]
            mask.dens[i, ] <- dnorm(ss.capt[i], mu.ss, det.pars[["sigma.ss"]])
        } else {
            stop("The binary capture history must only contain 0s and 1s.")
        }
    }
    apply(mask.dens, 2, prod)
}

toa.density <- function(bincapt, covariates, traps, mask, pars, dists, sound.speed){
  if (!('toa' %in% names(covariates))) {
    stop("Something went wrong calculating toa density,
         'toa' column missing from covariate matrix.")
  }
  
  dists <- dists[bincapt == 1, ]
  toa.capt <- covariates$toa[bincapt == 1]
  sigma.toa <- pars$sigma.toa
  prod.times <- toa.capt - dists / sound.speed
  toa.ssq <- apply(prod.times, 2, function(x) sum((x - mean(x))^2))
  out <- (2*pi*sigma.toa^2)^((1 - sum(bincapt))/2)*
    exp(toa.ssq/(-2*sigma.toa^2))
}

## Plots arrows on traps where a detection was made, showing estimated bearing.
show.arrows <- function(fit = NULL, id, session, arrow.length = NULL, trap.col, capt = NULL, traps = NULL){
    xlim <- par("usr")[c(1, 2)]
    ylim <- par("usr")[c(3, 4)]
    if (is.null(arrow.length)){
        arrow.length <- 0.05*min(c(diff(range(xlim)), diff(range(ylim))))
    }
    if (!is.null(fit)){
        capt.all <- get.capt(fit, session)
    } else {
        capt.all <- capt
    }
    capt <- capt.all$bincapt[id, ]
    bearing.capt <- capt.all$bearing[id, capt == 1]
    if (!is.null(fit)){
        trappos <- get.traps(fit, session)
    } else {
        trappos <- traps
    }
    trappos <- trappos[which(capt == 1), , drop = FALSE]
    sinb <- sin(bearing.capt)*arrow.length
    cosb <- cos(bearing.capt)*arrow.length
    arrows(trappos[, 1], trappos[, 2], trappos[, 1] + sinb, trappos[, 2] + cosb,
           length = 0.1, col = trap.col, lwd = 2)
}

## Plots circles around traps where a detection was made, showing estimated distance.
show.circles <- function(fit = NULL, id, session, trap.col, capt = NULL, traps = NULL){
    if (!is.null(fit)){
        capt.all <- get.capt(fit, session)
    } else {
        capt.all <- capt
    }
    capt <- capt.all$bincapt[id, ]
    dist.capt <- capt.all$dist[id, capt == 1]
    if (!is.null(fit)){
        trappos <- get.traps(fit, session)
    } else {
        trappos <- traps
    }
    trappos <- trappos[which(capt == 1), , drop = FALSE]
    for (i in 1:nrow(trappos)){
        centre <- trappos[i, ]
        radius <- dist.capt[i]
        circles(centre, radius, col = trap.col, lwd = 2)
    }
}

circles <- function(centre, radius, ...){
    bearings <- seq(0, 2*pi, length.out = 100)
    xs <- centre[1] + sin(bearings)*radius
    ys <- centre[2] + cos(bearings)*radius
    lines(xs, ys, ...)
}
