## Unsuitable habitat: mask points where there can be no animals (or calls), which are removed
## from the mask. The 'unsuitable' arguments of read.acre() and fit.acre() are logical
## expressions in the mask-level covariates and the mask coordinates, x and y, e.g.,
## unsuitable = (FOREST_COVER == "NOT_DENSE") | (VILLAGE < 1). Each is kept as a one-sided
## formula, so it can be evaluated again later (e.g., for a new mask, or to predict density) in
## the environment where it was written.

## The 'unsuitable' argument as a one-sided formula.
##
## expr: The unevaluated argument, from substitute(). It can also be a one-sided formula, or the
##       name of an object that holds one, for use in code.
## env:  The environment the argument was provided in.
unsuitable.formula <- function(expr, env){
  if(is.null(expr)) return(NULL)
  if(is.symbol(expr) && exists(as.character(expr), envir = env)){
    obj <- get(as.character(expr), envir = env)
    if(inherits(obj, 'formula')) expr <- obj
  }
  if(inherits(expr, 'formula') || (is.call(expr) && identical(expr[[1]], as.name('~')))){
    f <- eval(expr, env)
    if(length(f) != 2) stop("A formula for 'unsuitable' must be one-sided, e.g., ~ VILLAGE < 1.")
    return(f)
  }
  stats::as.formula(call('~', expr), env = env)
}

## Which points are in unsuitable habitat.
##
## unsuitable: A list of formulas from unsuitable.formula(). A point is unsuitable if any of
##             them is TRUE for it.
## data:       A data frame with columns x and y, and covariates, with a row for each point.
unsuitable.points <- function(unsuitable, data){
  out <- rep(FALSE, nrow(data))
  for(f in unsuitable){
    expr <- f[[2]]
    u <- tryCatch(eval(expr, data, environment(f)),
                  error = function(e) stop("Cannot evaluate 'unsuitable' (", deparse1(expr), "): ",
                                           conditionMessage(e), call. = FALSE))
    if(!is.logical(u) || !(length(u) %in% c(1, nrow(data)))){
      stop("'unsuitable' (", deparse1(expr), ") must give TRUE or FALSE for each mask point.", call. = FALSE)
    }
    if(anyNA(u)){
      stop("'unsuitable' (", deparse1(expr), ") is NA for some mask points.", call. = FALSE)
    }
    out <- out | u
  }
  out
}

## Removes unsuitable habitat from the masks, and from the mask-level covariates.
##
## mask:       A list of masks, one for each session.
## mask.cov:   A data frame of mask-level covariates, with columns 'session' and 'mask', or NULL.
## unsuitable: A list of formulas from unsuitable.formula().
##
## Returns a list with the new masks, and the new covariates, with mask points renumbered.
remove.unsuitable <- function(mask, mask.cov, unsuitable){
  if(is.matrix(mask) || is.data.frame(mask)) mask <- list(mask)
  keep.cov <- rep(TRUE, NROW(mask.cov))
  for(s in seq_along(mask)){
    m <- mask[[s]]
    data <- data.frame(x = m[, 1], y = m[, 2])
    if(!is.null(mask.cov)){
      rows <- which(mask.cov$session == s)
      rows <- rows[order(mask.cov$mask[rows])]
      stopifnot(length(rows) == nrow(m))
      cov.names <- setdiff(colnames(mask.cov), c('session', 'mask', 'x', 'y'))
      data <- cbind(data, mask.cov[rows, cov.names, drop = FALSE])
    }
    bad <- unsuitable.points(unsuitable, data)
    if(all(bad)) stop("All mask points for session ", s, " are in unsuitable habitat.")
    area <- attr(m, 'area')
    buffer <- attr(m, 'buffer')
    m <- m[!bad, , drop = FALSE]
    attr(m, 'area') <- area
    attr(m, 'buffer') <- buffer
    mask[[s]] <- m
    if(!is.null(mask.cov)){
      keep.cov[rows[bad]] <- FALSE
      mask.cov$mask[rows[!bad]] <- seq_len(sum(!bad))
    }
  }
  if(!is.null(mask.cov)){
    mask.cov <- mask.cov[keep.cov, , drop = FALSE]
    rownames(mask.cov) <- NULL
  }
  list(mask = mask, mask.cov = mask.cov)
}
