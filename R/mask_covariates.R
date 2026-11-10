## Covariates attached to masks: a mask provided to read.acre() or fit.acre() can be a data frame
## (or a list of them, one for each session) with columns x and y for the mask points, and other
## columns for covariate values at those points. These are used as they are, rather than being
## interpolated, and in favour of any covariates with the same names provided elsewhere.

## Splits covariates attached to masks from the mask coordinates.
##
## mask: A list of masks, one for each session, each a matrix or data frame.
##
## Returns a list with the masks, with only their coordinates (and their attributes), and a data
## frame of the covariates, with columns session and mask, or NULL if there are none.
split.mask.cov <- function(mask){
  cov_names <- lapply(mask, mask.cov.names)
  if(all(lengths(cov_names) == 0)) return(list(mask = lapply(mask, mask.coords), cov = NULL))
  if(any(vapply(cov_names, function(x) !setequal(x, cov_names[[1]]), logical(1)))){
    stop("Masks must have the same covariates attached for every session.", call. = FALSE)
  }
  cov_names <- cov_names[[1]]
  if(any(c("session", "mask") %in% cov_names)){
    stop("Covariates attached to a mask cannot be called 'session' or 'mask'.", call. = FALSE)
  }
  cov <- vector("list", length(mask))
  for(s in seq_along(mask)){
    m <- mask[[s]]
    df <- as.data.frame(m, stringsAsFactors = FALSE)
    mask[[s]] <- mask.coords(m)
    cov[[s]] <- data.frame(session = s, mask = seq_len(nrow(df)), df[, cov_names, drop = FALSE],
                           stringsAsFactors = FALSE)
    ## Categorical covariates are characters, as they are when interpolated.
    for(i in cov_names) if(is.factor(cov[[s]][[i]])) cov[[s]][[i]] <- as.character(cov[[s]][[i]])
    ## Missing values are not interpolated; that is what loc.cov is for.
    n_na <- vapply(cov[[s]][cov_names], function(x) sum(is.na(x)), numeric(1))
    if(any(n_na > 0)){
      stop("Covariates attached to a mask cannot have missing values, but ",
           paste0("'", names(n_na)[n_na > 0], "' has ", n_na[n_na > 0], collapse = ", "),
           " for session ", s, ". Use 'loc.cov' to interpolate covariates at mask points instead.",
           call. = FALSE)
    }
  }
  cov <- do.call(rbind, cov)
  rownames(cov) <- NULL
  list(mask = mask, cov = cov)
}

## The coordinates of a mask's points, as a matrix with columns x and y, keeping its area and
## buffer attributes. A mask with only two columns can have them unnamed.
mask.coords <- function(m){
  coords <- if(NCOL(m) > 2) as.matrix(as.data.frame(m)[, c("x", "y")]) else as.matrix(m)
  storage.mode(coords) <- "double"
  colnames(coords) <- c("x", "y")
  for(a in c("area", "buffer")) attr(coords, a) <- attr(m, a)
  coords
}

## The names of the covariates attached to a mask: its columns other than x and y.
##
## mask: A mask, as a matrix or data frame, or a list of them, one for each session.
mask.cov.names <- function(mask){
  if(is.list(mask) && !is.data.frame(mask)) mask <- mask[[1]]
  if(is.null(mask) || NCOL(mask) <= 2) return(character(0))
  if(!all(c("x", "y") %in% colnames(mask))){
    stop("A mask with covariates attached must have columns 'x' and 'y' for the mask points.", call. = FALSE)
  }
  setdiff(colnames(mask), c("x", "y"))
}

## Removes covariates from loc.cov, dist.cov, time.loc.cov, session.cov and trap.cov that have the
## same names as covariates attached to a mask, with a warning.
##
## names: The names of the covariates attached to the mask.
## covs:  A list with components loc.cov, dist.cov, time.loc.cov, session.cov and trap.cov, any
##        of which can be NULL.
##
## Returns covs, with these covariates removed. Components with no covariates left are NULL.
drop.mask.cov.clashes <- function(names, covs){
  warn <- function(found, source){
    if(length(found) > 0){
      warning("Covariate", if(length(found) > 1) "s", " ", paste0("'", found, "'", collapse = ", "),
              " in '", source, "' ", if(length(found) > 1) "are" else "is",
              " ignored in favour of the covariate", if(length(found) > 1) "s",
              " attached to the mask.", call. = FALSE)
    }
  }
  ## Data frames, or lists of them, with coordinates (and times) as well as covariates.
  for(source in c("loc.cov", "time.loc.cov")){
    x <- covs[[source]]
    if(is.null(x)) next
    keep <- c("x", "y", "time")
    lst <- if(is.data.frame(x)) list(x) else x
    warn(intersect(unique(unlist(lapply(lst, colnames))), names), source)
    lst <- lapply(lst, function(df) df[, setdiff(colnames(df), names), drop = FALSE])
    lst <- Filter(function(df) length(setdiff(colnames(df), keep)) > 0, lst)
    covs[source] <- list(if(length(lst) == 0) NULL else if(is.data.frame(x)) lst[[1]] else lst)
  }
  ## Features named after the covariates.
  if(!is.null(covs$dist.cov)){
    warn(intersect(names(covs$dist.cov), names), "dist.cov")
    dist.cov <- covs$dist.cov[setdiff(names(covs$dist.cov), names)]
    covs["dist.cov"] <- list(if(length(dist.cov) == 0) NULL else dist.cov)
  }
  ## Session- and detector-level covariates.
  for(source in c("session.cov", "trap.cov")){
    x <- covs[[source]]
    if(is.null(x)) next
    found <- setdiff(intersect(colnames(x), names), c("session", "trap"))
    warn(found, source)
    x <- x[, setdiff(colnames(x), found), drop = FALSE]
    covs[source] <- list(if(length(setdiff(colnames(x), c("session", "trap"))) == 0) NULL else x)
  }
  covs
}

## Combines interpolated mask-level covariates with those attached to the mask. Interpolated
## covariates with the same names as attached ones are ignored, with a warning.
##
## interp:   A data frame of interpolated covariates from location.cov.to.mask(), or NULL.
## attached: A data frame of attached covariates from split.mask.cov(), or NULL.
combine.mask.cov <- function(interp, attached){
  if(is.null(attached)) return(interp)
  if(is.null(interp)) return(attached)
  att_names <- setdiff(colnames(attached), c("session", "mask"))
  clash <- intersect(att_names, colnames(interp))
  if(length(clash) > 0){
    drop.mask.cov.clashes(att_names, list(loc.cov = interp[, c("session", "mask", clash)]))
    interp <- interp[, setdiff(colnames(interp), clash), drop = FALSE]
  }
  out <- merge(interp, attached, by = c("session", "mask"), sort = FALSE)
  out <- out[order(out$session, out$mask), , drop = FALSE]
  rownames(out) <- NULL
  out
}
