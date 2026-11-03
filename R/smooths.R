## Smooth terms in the linear predictor of an extended parameter.
##
## The design matrix of a smooth term depends on the data it was built from (e.g., the basis of a
## thin plate regression spline), so design matrices for new data must be built from the smooth
## objects created when the model was fitted, using mgcv::PredictMat(), rather than by setting up
## the model formula again on the new data.
##
## A penalised smooth is represented as a mixed model (Wood, 2017, Section 5.8; Seaton et al., in
## submission): mgcv::smooth2random() reparameterises its basis so that it splits into columns for
## the unpenalised null space of the penalty, whose coefficients are fixed effects, and columns
## for the penalised range space, whose coefficients are independent normal random effects with
## a common standard deviation. TMB integrates the random effects out of the likelihood using the
## Laplace approximation, and the standard deviation controls how wiggly the smooth is.

## Information required to build the design matrix for a mask-level linear predictor from any
## data, with the columns for fixed coefficients first, then the columns for the random effects.
##
## G:        A gam() setup object (from mgcv::gam(fit = FALSE)) for the mask-level formula.
## data:     The data used to set up G.
## penalise: Whether or not the smooth terms are penalised. Smooths with fx = TRUE are never
##           penalised.
smooth.info <- function(G, data, penalise){
  ## The parametric part of the model, which includes the intercept.
  pf <- mgcv::interpret.gam(G$formula)$pf
  G_param <- mgcv::gam(pf, data = data, fit = FALSE)
  smooths <- vector('list', length(G$smooth))
  for(j in seq_along(G$smooth)){
    sm <- G$smooth[[j]]
    re <- NULL
    if(penalise && !sm$fixed){
      sm$X <- G$X[, sm$first.para:sm$last.para, drop = FALSE]
      re <- tryCatch(mgcv::smooth2random(sm, "", type = 2),
                     error = function(e) stop("The smooth term ", sm$label, " cannot be penalised: ",
                                              conditionMessage(e), call. = FALSE))
      sm$X <- NULL
    }
    smooths[[j]] <- list(sm = sm, re = re)
  }
  list(G_param = G_param, smooths = smooths)
}

## The design matrix for a mask-level linear predictor, built from a smooth.info() object.
##
## info:    An object returned by smooth.info().
## newdata: A data frame of covariates.
## X.param: The design matrix for the parametric part of the model, if it is already available
##          (i.e., if newdata are the data used to set up the model).
##
## Returns a list with components X, the columns for fixed coefficients (including the intercept),
## Z, the columns for random effects, and re_fac, which gives, for each column of Z, the index of
## the standard deviation of its random effect. The standard deviations are labelled by sd_labels.
smooth.design <- function(info, newdata, X.param = NULL){
  X <- if(is.null(X.param)) param.design.new(info$G_param, newdata) else X.param
  colnames(X) <- colnames(info$G_param$X)
  Z <- matrix(0, nrow = nrow(X), ncol = 0)
  re_fac <- integer(0)
  sd_labels <- character(0)
  for(s in info$smooths){
    cols <- smooth.columns(s, newdata)
    X <- cbind(X, cols$X)
    Z <- cbind(Z, cols$Z)
    re_fac <- c(re_fac, length(sd_labels) + cols$re_fac)
    sd_labels <- c(sd_labels, cols$sd_labels)
  }
  list(X = X, Z = Z, re_fac = re_fac, sd_labels = sd_labels)
}

## The design matrix columns for one smooth, for new data, as in smooth.design(), but with re_fac
## indexing this smooth's standard deviations only.
##
## s:       A component of the smooths element of a smooth.info() object.
## newdata: A data frame of covariates.
smooth.columns <- function(s, newdata){
  sm <- s$sm
  if(is.null(s$re)){
    X <- mgcv::PredictMat(sm, newdata)
    colnames(X) <- paste(sm$label, seq_len(ncol(X)), sep = ".")
    return(list(X = X, Z = matrix(0, nrow = nrow(X), ncol = 0), re_fac = integer(0),
                sd_labels = character(0)))
  }
  r <- smooth2random.predict(sm, s$re, newdata)
  X <- r$Xf
  if(ncol(X) > 0) colnames(X) <- paste0(sm$label, ".f", seq_len(ncol(X)))
  Z <- matrix(0, nrow = nrow(X), ncol = 0)
  re_fac <- integer(0)
  sd_labels <- character(0)
  n_rand <- length(r$rand)
  for(b in seq_len(n_rand)){
    Z_b <- r$rand[[b]]
    sd_label <- if(n_rand == 1) sm$label else paste(sm$label, b, sep = ".")
    colnames(Z_b) <- paste0(sd_label, ".r", seq_len(ncol(Z_b)))
    Z <- cbind(Z, Z_b)
    sd_labels <- c(sd_labels, sd_label)
    re_fac <- c(re_fac, rep(length(sd_labels), ncol(Z_b)))
  }
  list(X = X, Z = Z, re_fac = re_fac, sd_labels = sd_labels)
}

## The terms in the model for D whose effect is plotted by plot(fit, type = "effect").
##
## fit:  A fitted model object.
## term: A smooth's label (e.g., "s(VILLAGE)"), for that smooth alone, or the name of a covariate
##       (e.g., "VILLAGE"), for all the terms that involve it: smooths, and parametric terms such
##       as VILLAGE and I(VILLAGE^2).
##
## Returns a list with the selected smooths (components of a smooth.info() object), the columns of
## the parametric design matrix for the selected parametric terms, the covariates the effect
## depends on, and a label.
effect.terms <- function(fit, term){
  info <- if('D' %in% get.par.extend.name(fit)) get.gam(fit, 'D')$gam_mask$smooth_info
  if(is.null(info)){
    stop("Effects can only be plotted for models for D with mask-level covariates.", call. = FALSE)
  }
  smooths <- info$smooths
  names(smooths) <- vapply(smooths, function(s) s$sm$label, character(1))
  G <- info$G_param
  p_labels <- attr(G$pterms, "term.labels")
  p_vars <- lapply(p_labels, function(l) all.vars(str2lang(l)))
  if(term %in% names(smooths)){
    sel_s <- smooths[term]
    sel_p <- integer(0)
  } else {
    sel_s <- Filter(function(s) term %in% s$sm$term, smooths)
    sel_p <- which(vapply(p_vars, function(v) term %in% v, logical(1)))
  }
  if(length(sel_s) + length(sel_p) == 0){
    stop("'", term, "' is not a term or covariate in the model for D.", call. = FALSE)
  }
  if(any(lengths(p_vars[sel_p]) > 1)){
    stop("Effects of covariates in interactions (", paste(p_labels[sel_p][lengths(p_vars[sel_p]) > 1], collapse = ", "),
         ") are not available.", call. = FALSE)
  }
  if(any(!vapply(sel_s, function(s) identical(s$sm$by, "NA"), logical(1)))){
    stop("Effects of smooths with a 'by' variable are not available.", call. = FALSE)
  }
  vars <- unique(c(unlist(p_vars[sel_p]), unlist(lapply(sel_s, function(s) s$sm$term))))
  if(any(vapply(G$var.summary[intersect(vars, names(G$var.summary))], is.factor, logical(1)))){
    stop("Effects of factors are not available.", call. = FALSE)
  }
  if(length(vars) > 2){
    stop("Effects of more than two covariates are not available.", call. = FALSE)
  }
  list(smooths = sel_s, cols = colnames(G$X)[G$assign %in% sel_p], vars = vars,
       label = paste(c(p_labels[sel_p], names(sel_s)), collapse = " + "), G_param = G)
}

## The covariates at the mask points used to fit the model, for all sessions, including the mask
## coordinates, in the order of the masks.
effect.mask.data <- function(fit){
  mask_data <- get.par.extend.data(fit)$mask
  mask_data <- mask_data[order(mask_data$session, mask_data$mask), , drop = FALSE]
  coords <- do.call(rbind, lapply(get.mask(fit), function(m) m[, 1:2, drop = FALSE]))
  for(i in setdiff(c("x", "y"), colnames(mask_data))) mask_data[[i]] <- coords[, i]
  mask_data
}

## The estimated effect of terms in the model for D, on the scale of the linear predictor, with
## standard errors from the covariance of their coefficients (including random effects) only.
## Smooths are centred by their identifiability constraints, like mgcv's plot.gam() shows them,
## so they average zero over the mask points used to fit the model, and parametric terms are
## centred to do the same, like termplot().
##
## fit:     A fitted model object.
## term:    As for effect.terms().
## newdata: A data frame with the covariates the effect depends on, or, if it depends on one
##          covariate, a vector of its values.
##
## Returns a data frame with those covariates, and columns est and se.
term.effect <- function(fit, term, newdata){
  et <- effect.terms(fit, term)
  if(!is.data.frame(newdata)){
    if(length(et$vars) > 1) stop("'newdata' must be a data frame with ", paste(et$vars, collapse = " and "), ".")
    newdata <- stats::setNames(data.frame(newdata), et$vars)
  }
  if(!all(et$vars %in% names(newdata))){
    stop("'newdata' must include ", paste(setdiff(et$vars, names(newdata)), collapse = " and "), ".", call. = FALSE)
  }
  X <- matrix(0, nrow = nrow(newdata), ncol = 0)
  if(length(et$cols) > 0){
    ## The other parametric covariates are held at their values at the first mask point, which do
    ## not affect the selected columns.
    nd <- effect.mask.data(fit)[rep(1, nrow(newdata)), , drop = FALSE]
    nd[et$vars] <- newdata[et$vars]
    X_p <- param.design.new(et$G_param, nd)
    colnames(X_p) <- colnames(et$G_param$X)
    X_p <- sweep(X_p[, et$cols, drop = FALSE], 2, colMeans(et$G_param$X[, et$cols, drop = FALSE]))
    X <- cbind(X, X_p)
  }
  for(s in et$smooths){
    cols <- smooth.columns(s, newdata)
    X <- cbind(X, cols$X, cols$Z)
  }
  coef_names <- paste0("D.", colnames(X), "_link")
  b <- coef(fit, types = "linked", pars = "D", re = TRUE)
  b <- stats::setNames(as.vector(b), if(is.null(dim(b))) names(b) else rownames(b))
  V <- vcov(fit, types = "linked", pars = "D", re = TRUE)
  if(is.list(V)) V <- V[[1]]
  V <- V[coef_names, coef_names, drop = FALSE]
  out <- newdata[, et$vars, drop = FALSE]
  rownames(out) <- NULL
  out$est <- drop(X %*% b[coef_names])
  out$se <- sqrt(pmax(rowSums((X %*% V) * X), 0))
  out
}

## Plots the estimated effects of terms in the model for D, on the scale of the linear predictor,
## from term.effect(). An effect of one covariate is a curve, with a shaded band (or dashed lines)
## for a confidence interval. An effect of two covariates (e.g., s(x, y)) is a surface over space,
## optionally followed by a surface of its standard errors. Returns the plotted effects invisibly,
## as a list of data frames from term.effect(), named by term.
##
## term:      Smooth labels or covariate names, as for effect.terms(); by default, all smooths.
## newdata:  For effects of two covariates, a data frame of points at which to show it, with x,
##            y, and the covariates. By default, the pixels of an x.pixels by y.pixels grid that
##            are within the masks' buffers of the detectors, including unsuitable habitat (or,
##            with p.dot.min, that contain detectable mask points), which is only available if
##            the covariates are x and y.
## n:         The number of covariate values at which to estimate an effect of one covariate.
## se:        Whether to show a confidence interval for an effect of one covariate (TRUE by
##            default), or a surface of standard errors for one of two (FALSE by default).
## shade:     Whether to show a confidence interval as a shaded band, rather than dashed lines.
## level:     The confidence level.
## p.dot.min: If provided, only mask points with a detection probability above this, from p.dot(),
##            set the range of a covariate, or the grid for a surface.
## add:       Whether to add a curve, or contours of a surface, to an existing plot.
## ask:       Whether to ask before showing each plot after the first.
## ...:       Other arguments for plot() or fields::imagePlot(), e.g., xlim, ylim, main, and for
##            curves, col, lty and lwd.
show.effect <- function(fit, term = NULL, newdata = NULL, n = 100, se = NULL, shade = TRUE,
                        level = 0.95, p.dot.min = NULL, add = FALSE, x.pixels = 100,
                        y.pixels = 100, ask = TRUE, ...){
  if(is.null(term)){
    info <- if('D' %in% get.par.extend.name(fit)) get.gam(fit, 'D')$gam_mask$smooth_info
    term <- vapply(info$smooths, function(s) s$sm$label, character(1))
    if(length(term) == 0) stop("The model for D has no smooths; use 'term' to select a covariate.", call. = FALSE)
  }
  ## The mask points that set the ranges of covariates.
  mask_data <- effect.mask.data(fit)
  if(!is.null(p.dot.min)){
    detectable <- unlist(p.dot(fit)) > p.dot.min
    if(!any(detectable)){
      stop("No mask points have a detection probability above p.dot.min = ", p.dot.min, ".", call. = FALSE)
    }
    mask_data <- mask_data[detectable, , drop = FALSE]
  }
  z_mult <- stats::qnorm(1 - (1 - level)/2)
  output <- list()
  n_plots <- 0
  for(j in seq_along(term)){
    et <- effect.terms(fit, term[j])
    plot.args <- list(...)
    if(length(et$vars) == 1){
      ## A curve.
      covariate <- et$vars
      xlim <- if(is.null(plot.args$xlim)) range(mask_data[[covariate]]) else plot.args$xlim
      eff <- term.effect(fit, term[j], seq(xlim[1], xlim[2], length.out = n))
      x <- eff[[covariate]]
      show.se <- if(is.null(se)) TRUE else se
      lower <- eff$est - z_mult*eff$se
      upper <- eff$est + z_mult*eff$se
      col <- if(is.null(plot.args$col)) 1 else plot.args$col
      lty <- if(is.null(plot.args$lty)) 1 else plot.args$lty
      lwd <- if(is.null(plot.args$lwd)) 1 else plot.args$lwd
      plot.args[c("col", "lty", "lwd")] <- NULL
      if(!add){
        if(ask && n_plots > 0) ask.user.for.next.plot()
        n_plots <- n_plots + 1
        if(is.null(plot.args$xlab)) plot.args$xlab <- covariate
        if(is.null(plot.args$ylab)) plot.args$ylab <- et$label
        if(is.null(plot.args$ylim)) plot.args$ylim <- range(if(show.se) c(lower, upper) else eff$est)
        do.call(graphics::plot, c(list(x = x, y = eff$est, type = "n"), plot.args))
      }
      if(show.se){
        if(shade){
          graphics::polygon(c(x, rev(x)), c(upper, rev(lower)), border = NA,
                            col = grDevices::adjustcolor(col, alpha.f = 0.25))
        } else {
          graphics::lines(x, upper, lty = 2, col = col, lwd = lwd)
          graphics::lines(x, lower, lty = 2, col = col, lwd = lwd)
        }
      }
      graphics::lines(x, eff$est, col = col, lty = lty, lwd = lwd)
    } else {
      ## A surface over space.
      if(!is.null(newdata)){
        points <- as.data.frame(newdata)
        if(!all(c("x", "y") %in% names(points))) stop("'newdata' must include x and y.", call. = FALSE)
      } else {
        if(!all(et$vars %in% c("x", "y"))){
          stop("Provide 'newdata' with ", paste(et$vars, collapse = " and "), " to show the effect of ",
               et$label, " over space.", call. = FALSE)
        }
        traps <- lapply(get.trap(fit), as.matrix)
        buffers <- lapply(get.mask(fit), attr, "buffer")
        if(is.null(p.dot.min) && !any(vapply(buffers, is.null, logical(1)))){
          ## Pixels within the masks' buffers of the detectors. These include unsuitable
          ## habitat, which is not in the masks, but where the effect is still defined.
          xlim <- range(unlist(mapply(function(tr, b) range(tr[, 1]) + c(-b, b), traps, buffers)))
          ylim <- range(unlist(mapply(function(tr, b) range(tr[, 2]) + c(-b, b), traps, buffers)))
          grid <- expand.grid(x = seq(xlim[1], xlim[2], length.out = x.pixels),
                              y = seq(ylim[1], ylim[2], length.out = y.pixels))
          keep <- rep(FALSE, nrow(grid))
          for(s in seq_along(traps)){
            d2 <- outer(grid$x, traps[[s]][, 1], "-")^2 + outer(grid$y, traps[[s]][, 2], "-")^2
            keep <- keep | apply(d2, 1, min) <= buffers[[s]]^2
          }
          points <- grid[keep, , drop = FALSE]
        } else {
          ## Pixels that contain (detectable) mask points.
          gx <- seq(min(mask_data$x), max(mask_data$x), length.out = x.pixels)
          gy <- seq(min(mask_data$y), max(mask_data$y), length.out = y.pixels)
          ix <- findInterval(mask_data$x, c(-Inf, (gx[-1] + gx[-x.pixels])/2, Inf))
          iy <- findInterval(mask_data$y, c(-Inf, (gy[-1] + gy[-y.pixels])/2, Inf))
          cells <- unique(cbind(ix, iy))
          points <- data.frame(x = gx[cells[, 1]], y = gy[cells[, 2]])
        }
      }
      eff <- term.effect(fit, term[j], points)
      if(!all(c("x", "y") %in% names(eff))) eff <- cbind(points[, c("x", "y")], eff)
      ux <- sort(unique(eff$x))
      uy <- sort(unique(eff$y))
      show.se <- if(is.null(se)) FALSE else se
      surfaces <- if(show.se && !add) c("est", "se") else "est"
      for(k in surfaces){
        z <- matrix(NA, nrow = length(ux), ncol = length(uy))
        z[cbind(match(eff$x, ux), match(eff$y, uy))] <- eff[[k]]
        if(add){
          graphics::contour(ux, uy, z, add = TRUE, ...)
          next
        }
        if(ask && n_plots > 0) ask.user.for.next.plot()
        n_plots <- n_plots + 1
        args <- plot.args
        if(k == "est"){
          zmax <- max(abs(eff$est))
          if(is.null(args$zlim)) args$zlim <- c(-zmax, zmax)
          if(is.null(args$col)) args$col <- grDevices::hcl.colors(100, "Blue-Red 3")
          if(is.null(args$main)) args$main <- et$label
        } else {
          args$zlim <- NULL
          args$col <- grDevices::hcl.colors(100, "viridis")
          args$main <- paste("Standard error:", et$label)
        }
        if(is.null(args$xlab)) args$xlab <- "x"
        if(is.null(args$ylab)) args$ylab <- "y"
        if(is.null(args$asp)) args$asp <- 1
        do.call(fields::imagePlot, c(list(x = ux, y = uy, z = z), args))
      }
    }
    output[[term[j]]] <- eff
  }
  invisible(output)
}

## The design matrix for the parametric part of a model, for new data. Parametric terms do not
## depend on the data used to set up the model, except for the levels of factors, which
## get.DX.new.gam() handles.
param.design.new <- function(G_param, newdata){
  if(length(G_param$var.summary) == 0){
    matrix(1, nrow = nrow(newdata), ncol = 1)
  } else {
    get.DX.new.gam(G_param, newdata)
  }
}

## Design matrices for new data under the random-effects reparameterisation of a smooth, from
## mgcv's documentation for smooth2random().
##
## sm:   A smooth object.
## re:   The result of mgcv::smooth2random(sm, ..., type = 2).
## data: A data frame of covariates.
smooth2random.predict <- function(sm, re, data){
  X <- mgcv::PredictMat(sm, data)
  ## Transform to the random-effects parameterisation.
  if(!is.null(re$trans.U)) X <- X %*% re$trans.U
  X <- t(t(X) * re$trans.D)
  ## Reorder the columns according to the random-effects reordering.
  X[, re$rind] <- X[, re$pen.ind != 0]
  ## Reorder the penalty index in the same way.
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for(i in seq_along(re$rand)){
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
  }
  r
}

## The design matrix for the mask-level part of an extended parameter's linear predictor, for new
## data: the columns of smooth.design() without the intercept, in the order used to fit the model.
mask.design.new <- function(gam.mask, newdata){
  if(is.null(gam.mask$smooth_info)){
    ## Fits from earlier versions of the package.
    return(get.DX.new.gam(gam.mask, newdata)[, -1, drop = FALSE])
  }
  d <- smooth.design(gam.mask$smooth_info, newdata)
  cbind(d$X[, -1, drop = FALSE], d$Z)
}

## Smooths are penalised unless they have fx = TRUE, but penalisation is only available for
## mask-level covariates in the model for D, so stops if any other smooth is not fixed.
##
## G:     A gam() setup object.
## par:   The name of the parameter.
## where: A description of the smooths that cannot be penalised, for the error message.
check.smooths.fixed <- function(G, par, where){
  penalised <- !vapply(G$smooth, function(sm) isTRUE(sm$fixed), logical(1))
  if(any(penalised)){
    labels <- vapply(G$smooth[penalised], function(sm) sm$label, character(1))
    stop("Penalised smooths are not available for ", where, " (", paste(labels, collapse = ", "),
         " in the model for ", par, "); use s(..., fx = TRUE) for an unpenalised smooth.", call. = FALSE)
  }
}

## Removes the predicted random effects for penalised smooths from the output of coef(), vcov(),
## stdEr() or confint(). They are not parameters, so they are reported in the fit's D.re component
## instead. Internal callers that need them alongside D's coefficients (e.g., to predict density,
## or to record bootstrap estimates) pass re = TRUE, which reaches this function through "...".
drop.re <- function(x, object, re = FALSE, ...){
  if(isTRUE(re) || is.null(object$D.smooth)) return(x)
  re_names <- paste0("D.", object$output.tmb$D_re_names)
  re_names <- c(re_names, paste0(re_names, "_link"))
  cl <- oldClass(x)
  if(is.list(x)){
    x[] <- lapply(x, drop.re, object = object)
  } else if(is.matrix(x)){
    keep <- !rownames(x) %in% re_names
    ## Covariance matrices lose columns as well as rows.
    if(identical(rownames(x), colnames(x))){
      x <- x[keep, keep, drop = FALSE]
    } else {
      x <- x[keep, , drop = FALSE]
    }
  } else {
    x <- x[!names(x) %in% re_names]
  }
  oldClass(x) <- cl
  x
}
