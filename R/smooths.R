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
    sm <- s$sm
    if(is.null(s$re)){
      X_sm <- mgcv::PredictMat(sm, newdata)
      colnames(X_sm) <- paste(sm$label, seq_len(ncol(X_sm)), sep = ".")
      X <- cbind(X, X_sm)
    } else {
      r <- smooth2random.predict(sm, s$re, newdata)
      if(ncol(r$Xf) > 0){
        colnames(r$Xf) <- paste0(sm$label, ".f", seq_len(ncol(r$Xf)))
        X <- cbind(X, r$Xf)
      }
      n_rand <- length(r$rand)
      for(b in seq_len(n_rand)){
        Z_b <- r$rand[[b]]
        sd_label <- if(n_rand == 1) sm$label else paste(sm$label, b, sep = ".")
        colnames(Z_b) <- paste0(sd_label, ".r", seq_len(ncol(Z_b)))
        Z <- cbind(Z, Z_b)
        sd_labels <- c(sd_labels, sd_label)
        re_fac <- c(re_fac, rep(length(sd_labels), ncol(Z_b)))
      }
    }
  }
  list(X = X, Z = Z, re_fac = re_fac, sd_labels = sd_labels)
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
