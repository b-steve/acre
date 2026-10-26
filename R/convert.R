#' Create mask object
#'
#' Creates a mask object to use with the function read.acre().
#'
#' @param traps a matrix or a data frame with two columns or a list of such matrices or data frames for a multi-session model.
#'              Each row in a matrix/data frame provides Cartesian coordinates (in metres) for the location of a detector.
#'              In a list of matrices or data frames, each element of the list corresponds to the detector location of a different session.
#'              If the detector locations stayed the same across several sessions, only one matrix/data frame is required.
#' @param buffer a scalar, the furthest distance that a detector could detect (in metres)
#' @param ...
#'
#' @export
create.mask <- function(traps, buffer, ...){
    ## Changing data frame to matrix to avoid is.list() issues.
    if (is(traps, 'data.frame')){
        traps <- as.matrix(traps)
    }
    if (is.list(traps)){
        ## Calling create.mask() on each session separately.
        n.sessions <- length(traps)
        mask <- vector(mode = "list", length = n.sessions)
        for (i in 1:n.sessions){
            mask[[i]] <- create.mask(traps[[i]], buffer = buffer, ...)
        }
    } else {
        ## Creating mask for a single session.
        traps <- convert.traps(traps)
        mask <- secr::make.mask(traps, buffer = buffer, type = "trapbuffer", ...)
        A <- attr(mask, "area")
        mask <- as.matrix(mask)
        attr(mask, "area") <- A
        attr(mask, "buffer") <- buffer
    }
    mask
}

#' Create captures object
#'
#' Create a captures object to use with the function read.acre
#'
#' @section Captures argument
#' The `captures` argument will be passed to the function `create.capt`. The main arguments provided are the session, the
#' identification number of the detected animal or sound, and the identification number of the trap which made the detection
#' (where the identification number is the row number of the corresponding trap in the matrix of trap locations. These columns
#' must be exactly labelled "session", "ID" and "trap" respectively.
#'
#' Additional columns can specify the auxiliary information collected over the course of the survey:
#'
#' \itemize{
#'  \item A column named `bearing`, containing estimated bearings (in radians) from the detector to each detected animal or sound.
#'
#'  \item A column named `dist`, containing the estimated distance between the detected animal or sound.
#'
#'  \item A column named `ss` containing the measured signal strengh of the detected sound.
#'
#'  \item A column named `toa` containing the measured time of arrival (in seconds) since the start of the survey (or some other
#' reference time) of the detected sound (only possible when the detectors are microphones).
#'
#'  \item A column named `animal_ID` containing the identification number of animals can be provided if individuals could be distinguished by
#'  their acoustic detection.
#'  }
#'
#' @param captures a data frame with capture data. Columns named "session", "ID" and "trap" are required, with
#'                 each row being regarded as one detection. Extra information can be provided as columns "bearing", "dist"
#'                 ,"ss", "toa", "animal_ID". See 'Captures argument'.
#' @param traps a matrix or a data frame with two columns or a list of such matrices or data frames for a multi-session model.
#'              Each row in a matrix/data frame provides Cartesian coordinates (in metres) for the location of a detector.
#'              In a list of matrices or data frames, each element of the list corresponds to the detector location of a different session.
#'              If the detector locations stayed the same across all sessions, only one matrix/data frame is required.
#' @return A data frame of the detections, one row per detection, with the auxiliary data that were
#'         recorded; it is checked, and sorted by session, animal (if any), call and detector.
#' @keywords internal
create.capt <- function(captures, traps){
    ## Make sure captures is a data frame or matrix with at least the session, ID and trap columns,
    ## and traps is a list, matrix or data frame.
    stopifnot(is.data.frame(captures) | is.matrix(captures),
              ncol(captures) >= 3,
              any(is.list(traps), is.matrix(traps), is.data.frame(traps)))

    expected_col_names <- c('session', 'animal_ID', 'ID', 'trap', 'bincapt',
                            'bearing', 'toa', 'dist', 'ss', 'occasion')
    if (!all(names(captures) %in% expected_col_names)) {
      unexpected_columns <- which(!(names(captures) %in% expected_col_names))
      warning("Unexpected column '", paste(names(captures)[unexpected_columns], collapse="', '"), "' found in captures dataframe. Unexpected columns will be ignored.")
    }

    if(is.matrix(captures)) captures <- as.data.frame(captures, stringsAsFactors = FALSE)

    ## Make sure that "session" and "trap" are numeric.
    captures$session <- as.numeric(captures$session)
    stopifnot(all(captures$session > 0), all(captures$session %% 1 == 0))
    stopifnot(all(captures$trap > 0), all(captures$trap %% 1 == 0))
    captures <- captures[order(captures$session, captures$ID, captures$trap),]

    ##------------------------------------------------------------------------------------------------------
    ## Deal with n.sessions.

    ## Determine the maximum session number in "captures".
    tem.n.sessions.capt <- max(captures$session)

    ## A "traps" of length 1 is a special case: the number of sessions comes from "captures".
    if(is(traps, 'list') && length(traps) > 1){
        n.sessions <- length(traps)
    } else {
        n.sessions <- tem.n.sessions.capt
    }

    ## The maximum session number in the captures data should be <= n.sessions.
    if(tem.n.sessions.capt > n.sessions) {
        stop('The number of sessions based on "traps" does not match the captures data.')
    }

    if(!all(captures$session %in% seq(n.sessions))){
        stop('"session" must be assigned as successive positive integers begins from one.')
    }

    ##-----------------------------------------------------------------------------------------------
    ## Deal with n.traps.

    ## Generate n.traps based on the captures data frame.
    tem.n.traps.capt <- numeric(n.sessions)
    tem.n.traps.capt[unique(captures$session)] <- aggregate(captures$trap,
                                                           list(session = captures$session), max)$x

    ## Check if "trap" is labeled by natural numbers.
    for(i in 1:n.sessions){
        tem <- subset(captures, captures$session == i)
        if(nrow(tem) > 0){
            if(any(!tem$trap %in% 1:tem.n.traps.capt[i])) {
                stop('labels of traps must be successive natural numbers')
            }
        }
    }

    ## If 'traps' has length 1, the same detectors were used in every session.
    if(any(is(traps, 'list') & length(traps) == 1, is.matrix(traps), is.data.frame(traps))){
        tem.traps <- vector('list', n.sessions)
        if(is(traps, 'list')){
            for(i in 1:n.sessions) tem.traps[[i]] <- traps[[1]]
        } else {
            for(i in 1:n.sessions) tem.traps[[i]] <- traps
        }
        traps <- tem.traps
    }
    n.traps <- sapply(traps, nrow)

    ## The value of n.traps generated from the capture data should always be smaller than the one from "traps".
    if(any(tem.n.traps.capt > n.traps)){
        msg <- paste0('In session(s): (', paste(which(tem.n.traps.capt > n.traps), collapse = ","),
                     '), the number of traps based on "traps" does not match the captures data')
        stop(msg)
    }

    stopifnot(all(n.traps > 0))

    ##-----------------------------------------------------------------------------------------------------
    is.animalID <- "animal_ID" %in% colnames(captures)

    dupli <- duplicated(captures[,c("session", "animal_ID"[is.animalID], "ID", "trap")])
    if(any(dupli)){
        warning("Ignoring the duplicated observations")
    }
    captures <- captures[!dupli,]

    ## One row per detection, with the auxiliary data that were recorded.
    keep <- c("session", "animal_ID", "ID", "trap", "bearing", "dist", "ss", "toa")
    captures <- captures[, colnames(captures) %in% keep, drop = FALSE]
    captures <- data.sort(captures, 'data.full')
    rownames(captures) <- NULL
    captures
}

#' A helper function to obtain the distances of the nearest points from a data frame
#'
#' @param from a matrix or a data frame with columns "x" and "y" contains the coordinates of the start points.
#' @param to a matrix or a data frame with columns "x" and "y" contains the coordinates of the end points.
#' @param col.name a character, contains the new column name for the distances in the output, default is
#'                 dist_nearest.
#'
#' @return a data frame as the same as "from", but with an extra column with assigned column name. For each row,
#'         it contains the distance to the nearest point in the "to" data set.
#' @keywords internal
dist.nearest <- function(from, to, col.name = 'dist_nearest'){
  stopifnot(all(c('x', 'y') %in% colnames(from)))
  stopifnot(all(c('x', 'y') %in% colnames(to)))
  output <- as.data.frame(from)
  from <- from[, c('x', 'y')]
  to <- to[, c('x', 'y')]

  nearest.point <- RANN::nn2(to, from, k = 1)

  output[[col.name]] <- as.vector(nearest.point$nn.dists)
  return(output)
}

#' Convert traps object
#'
#' Converts an `acre` traps matrix to a `secr` traps
#' object.
#'
#' The returned object is suitable for use as the `traps`
#' argument of the function [make.capthist].
#'
#' @param ss Logical, set to `TRUE` if a signal strength
#'     detection function is to be used.
#' @param traps a matrix or a data frame, contains one session's detectors' coordinates
#'
#' @return An object of class `traps` comprising a data frame of
#'     x- and y-coordinates, the detector type ('single', 'multi',
#'     'proximity', 'count', 'polygon' etc.), and possibly other
#'     attributes.
#' @keywords internal
convert.traps <- function(traps, ss = FALSE){
    if (is(traps, "list")){
        stop("The convert.traps() function will only convert single-session trap objects.")
    }
    n.traps <- nrow(traps)
    colnames(traps) <- c("x", "y")
    traps.df <- data.frame(names = 1:n.traps, traps)
    detector <- ifelse(ss, "signal", "proximity")
    secr::read.traps(data = traps.df, detector = detector)
}

## Some sort of helper function. I don't know what it does.
convert.one.mask <- function(x, trap){
  stopifnot(any(is(x, 'data.frame'), is(x, 'matrix')))
  if(is(x, 'data.frame')){
    x <- as.matrix(x)
  }
  stopifnot(ncol(x) == 2)
  stopifnot(is.numeric(x))
  colnames(x) <- c('x', 'y')
  if(is.null(attr(x, "buffer"))){
    stopifnot(any(is(trap, 'data.frame'), is(trap, 'matrix')))
    trap <- as.matrix(trap)
    stopifnot(is.numeric(trap))
    colnames(trap) <- c('x', 'y')

    d <- distances(trap, x)
    d_closet_trap <- apply(d, 2, min)
    buffer <- max(d_closet_trap)
    attr(x, "buffer") <- buffer
  }
  if(is.null(attr(x, "area"))){
    sp <- as.matrix(dist(x))
    spacing <- apply(sp, 1, function(x) min(x[x > 0]))
    spacing <- median(spacing, na.rm = T)
    area <- spacing^2/10000
    attr(x, "area") <- area
  }
  return(x)
}

