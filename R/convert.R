#' Create mask object
#'
#' Creates a mask object to use with the function read.acre
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
#' @param ind.model a logical value indicates whether to include individual identification which should be recorded as "animal_ID",
#'                  a column in the argument "captures". By default, it will be NULL, so the model will determine it automatically
#'                  depending on whether this information is provided in "captures".
#' @param n.sessions a numeric value denotes the number of sessions. It will only be used when the argument "traps" is not provided,
#'                   or it is a list with one element, or it is a data frame, or it is a matrix.
#' @param n.traps a numeric vector denotes the number of traps in each session. It will only be used when the argument "traps" is not provided.
#'                If its length is 1, its value will be used for each session when we have multiple sessions. If its length is greater than 1,
#'                the length of it must match the number of sessions.
#' @param mrds.locs a list with length of n.sessions with data frames or matrices as components, each data.frame or matrix
#'                  contains two columns record the Cartesian coordinates of each call. If a session has no detection,
#'                  keep the corresponding component as NULL. If animal.model, then each component must be a data.frame
#'                  with 4 columns: "animal_ID", "ID", "mrds_x" and "mrds_y".
#'
#' @return an object ready for model fitting.
#' @export
#'
#' @examples
create.capt <- function(captures, traps = NULL, ind.model = NULL, n.sessions = NULL, n.traps = NULL,
                        mrds.locs = NULL){

    # Make sure captures is DF or matrix
    # Make sure required columns are included in captures (at least 4)
    # Make sure traps is either list, matrix, data frame or NULL
    # Make sure at least one of traps, or (n.traps & n.sessions) provided
    stopifnot(is.data.frame(captures) | is.matrix(captures),
              ncol(captures) >= 3,
              any(is.null(traps), is.list(traps), is.matrix(traps), is.data.frame(traps)),
              any(!is.null(traps), !is.null(n.traps) & !is.null(n.sessions)))

    if ((!missing(n.traps) | !missing(n.sessions)) & is.null(traps)){
        warning("Arguments 'n.traps' and 'n.sessions' are deprecated. Please provide the the 'traps' argument instead.
        Future versions of acre will require the 'traps' argument to be provided to the 'create.capt()' function.")
    }
  
    expected_col_names <- c('session', 'animal_ID', 'ID', 'trap', 'bincapt',
                            'bearing', 'toa', 'dist', 'ss', 'mrds_x', 'mrds_y', 'occasion')
    if (!all(names(captures) %in% expected_col_names)) {
      unexpected_columns <- which(!(names(captures) %in% expected_col_names))
      warning("Unexpected column '", paste(names(captures)[unexpected_columns], collapse="', '"), "' found in captures dataframe. Unexpected columns will be ignored.")
    }

    if(is.matrix(captures)) captures = as.data.frame(captures, stringsAsFactors = FALSE)

    #make it sure that "session" and "trap" are numeric
    captures$session = as.numeric(captures$session)
    stopifnot(all(captures$session > 0), all(captures$session %% 1 == 0))
    stopifnot(all(captures$trap > 0), all(captures$trap %% 1 == 0))
    captures = captures[order(captures$session, captures$ID, captures$trap),]

    #------------------------------------------------------------------------------------------------------
    #deal with n.sessions


    #determine the maximum "captures" session number, and assume this is # sessions we are using
    tem.n.sessions.capt = max(captures$session)

    #if we have traps
    if(!is.null(traps)){
        #generate n.sessions based on "traps" data
        tem.n.sessions.trap = ifelse(is(traps, 'list'), length(traps), 1)

        #length of traps is 1 is a special case
        if(tem.n.sessions.trap == 1){
            #if n.sessions is not provided, we assume it should be max(captures$session)
            if(is.null(n.sessions)) n.sessions = tem.n.sessions.capt
            #if n.sessions is not null, we just use it, no code needed
        } else {
            #if length of traps is not 1, we should use it as n.sessions
            if(!is.null(n.sessions)) warning("'traps' is provided, argument 'n.sessions' will be overwriten")
            n.sessions = tem.n.sessions.trap
        }
    }

    #if "traps" is null, the constraint as very beginning of this function has confirmed that
    #"n.sessions" will be provided, so just use it, no code needed

    #the maximum of session in captures data should be <= n.sessions
    if(tem.n.sessions.capt > n.sessions) {
        stop('The number of sessions based on "n.sessions" or "traps" can not match the captures data')
    }


    if(!all(captures$session %in% seq(n.sessions))){
        stop('"session" must be assigned as successive positive integers begins from one.')
    }


    #-----------------------------------------------------------------------------------------------

    #deal with n.traps

    #generate n.traps based on captures data frame
    tem.n.traps.capt = numeric(n.sessions)
    tem.n.traps.capt[unique(captures$session)] = aggregate(captures$trap,
                                                           list(session = captures$session), max)$x

    #check if "trap" is labeled by natural numbers
    for(i in 1:n.sessions){
        tem = subset(captures, captures$session == i)
        if(nrow(tem) > 0){
            if(any(!tem$trap %in% 1:tem.n.traps.capt[i])) {
                stop('labels of traps must be successive natural numbers')
            }
        }
    }


    if(!is.null(traps)){
        #if 'traps' is provided deal with the special case that its length is 1
        if(any(is(traps, 'list') & length(traps) == 1, is.matrix(traps), is.data.frame(traps))){
            tem.traps = vector('list', n.sessions)
            if(is(traps, 'list')){
                for(i in 1:n.sessions) tem.traps[[i]] = traps[[1]]
            } else {
                for(i in 1:n.sessions) tem.traps[[i]] = traps
            }
            traps = tem.traps
        }

        if(!is.null(n.traps)) warning("'traps' is provided, argument 'n.traps' will be overwriten")
        n.traps = sapply(traps, nrow)
    } else {
        #length(n.traps) == 1 is a special case
        if(length(n.traps) == 1){
            n.traps = rep(n.traps, n.sessions)
        } else {
            msg = paste0("The length of 'n.traps' should match the number of sessions: ", n.sessions)
            if(length(n.traps) != n.sessions) stop(msg)
        }
    }

    #n.traps generated from capture data should always smaller than it from "traps"
    if(any(tem.n.traps.capt > n.traps)){
        msg = paste0('In session(s): (', paste(which(tem.n.traps.capt > n.traps), collapse = ","),
                     '), the number of traps based on "n.traps" or "traps" cannot match the captures data')
        stop(msg)
    }

    stopifnot(all(n.traps > 0))

    #-----------------------------------------------------------------------------------------------------

    #deal with mrds.locs. "ID" will be dealt with together as they are closely related
    #here we require mrds.locs a list with length of n.sessions, if there is no detection for
    #any sessions, the user should leave that component to be NULL.
    is.animalID.included = "animal_ID" %in% colnames(captures)

    if(is.null(ind.model)){
      is.animalID = is.animalID.included
    } else {
      is.animalID = is.animalID.included & ind.model
    }

    #when animalID is included but the user assign not to fit individual embedded model,
    #exclude the animalID information from the capture history
    if(is.animalID.included & !is.animalID){
      #re-assign cue ID to make sure no duplication ID after removing animal_ID
      captures$ID = as.numeric(as.factor(paste(captures$session, captures$animal_ID, captures$ID, sep = "-")))
      captures[,-which(colnames(captures) == "animal_ID")]
    }

    is.mrds <- !is.null(mrds.locs)

    #initial check, for further check about each component of mrds.locs, will be processed later
    #because there could be NULL for some component of this list.

    if (is.mrds){
        if(is(mrds.locs, "list")) {
            if(n.sessions != length(mrds.locs)) {
                msg = paste0("mrds.locs must be a list with a length of n.sessions: ", n.sessions)
                stop(msg)
            }

            if(is.animalID){
                stopifnot(all(sapply(mrds.locs, function(x) is(x, 'data.frame') | is.null(x))))
            } else {
                stopifnot(all(sapply(mrds.locs, function(x) is(x, 'data.frame') | is(x, 'matrix') | is.null(x))))
            }

        } else {
            stopifnot(n.sessions == 1)

            if(is.animalID){
                stopifnot(is(mrds.locs, 'data.frame'))
            } else {
                stopifnot(any(is(mrds.locs, 'data.frame'), is(mrds.locs, 'matrix')))
            }

            #convert it to a list with one element only, to make it be in consistent format of multi-session case
            mrds.locs = list(mrds.locs)
        }

    }


    dupli = duplicated(captures[,c("session", "animal_ID"[is.animalID], "ID", "trap")])
    if(any(dupli)){
        warning("Ignoring the duplicated observations")
    }
    captures = captures[!dupli,]

    #if not animal_ID data, use previous version output, which is a list with 'bincapt' and etc.
    #and if animal_ID data, then output a data frame directly since there is no necessary for compatibility

    #since 'session' is already be guaranteed to be successive natural numbers, we could directly aggregate
    #captures based on 'session', the output's order will be correct.
    if(!is.animalID){
        #n_IDs will be used later to construct the matrices of output
        n_IDs = numeric(n.sessions)
        n_IDs[unique(captures$session)] = aggregate(captures$ID, list(session = captures$session),
                                                    function(x) length(unique(x)))$x
    } else {
        #"n_animals_call" only used later to check 'mrds.locs'
        n_animals_call = numeric(n.sessions)
        n_animals_call[unique(captures$session)] = aggregate(paste(captures$animal_ID, captures$ID, sep = "_"),
                                                             list(session = captures$session),
                                                             function(x) length(unique(x)))$x
        n_animals = numeric(n.sessions)
        n_animals[unique(captures$session)] = aggregate(captures$animal_ID, list(session = captures$session),
                                                        function(x) length(unique(x)))$x
    }

    #check mrds.locs with 'ID'/'animal_ID' together
    if(is.mrds){
        n_col_mrds = lapply(mrds.locs, ncol)
        n_col_mrds = lapply(n_col_mrds, function(x) ifelse(is.null(x), 0, x))
        n_col_mrds = do.call('c', n_col_mrds)

        n_row_mrds = lapply(mrds.locs, nrow)
        n_row_mrds = lapply(n_row_mrds, function(x) ifelse(is.null(x), 0, x))
        n_row_mrds = do.call('c', n_row_mrds)

        #check whether the number of 'animal_ID' or 'ID' for each session matches the number of rows in 'mrds.locs'
        #and check whether the 'animal_ID' or 'ID' are successive natural numbers
        if(is.animalID){

            stopifnot(all(n_row_mrds == n_animals_call))

            #merge 'mrds.locs' into captures, the Cartesian coordinates of mrds will be recorded as 'mrds_x' and 'mrds_y'
            #these column names matter, because they will be used in the model fitting function directly, so cannot be renamed
            for(s in 1:n.sessions){
                tem = mrds.locs[[s]]
                if(!is.null(tem)){
                    stopifnot(ncol(tem) == 4)
                    stopifnot(all(c('animal_ID', 'ID', 'mrds_x', 'mrds_y') %in% colnames(tem)))
                    tem$session = s
                }
                #if it is NULL, just keep it as NULL
            }
            mrds.locs = do.call('rbind', mrds.locs)
            captures = merge(captures, mrds.locs, by = c('session', 'animal_ID', 'ID'), all.x = TRUE)
            #captures = convert_natural_number(dat = captures, is.animalID = is.animalID, which.convert = "both")

        } else {
            stopifnot(all(n_row_mrds == n_IDs))
            is.natural_number = natural_number_check(captures$session, captures$ID)
            if(!is.natural_number) stop('As "mrds.locs" is provided, "ID" only accepts successive natural numbers.')
        }
    } #else {
#        if(is.animalID){
#            captures = convert_natural_number(dat = captures, is.animalID = is.animalID, which.convert = 'both')
#        } else {
#            captures = convert_natural_number(dat = captures, is.animalID = is.animalID, which.convert = 'ID')
#        }
#    }

    captures = sort.data(captures, 'data.full')

    #-----------------------------------------------------------------------------------------------------
    #all checks have been done, below we generate output list
    if(!is.animalID){
        all.types <- c("bearing", "dist", "ss", "toa")
        info.types <- all.types[all.types %in% colnames(captures)]
        out.list <- vector(mode = "list", length = n.sessions)

        for (s in 1:n.sessions){
            tem <- captures[captures$session == s,, drop = FALSE]

            id <- tem$ID
            n <- n_IDs[s]

            out <- vector(mode = "list", length = length(info.types) + 1 + as.numeric(is.mrds))
            #if animal_ID is used, we use data frame for 'bincapt', and we add two columns after
            #n.traps[s] columns, so that the order of the first n.traps[s] columns still make sense


            for (i in 1:length(out)){
                out[[i]] <- matrix(0, nrow = n, ncol = n.traps[s])
            }

            names(out) <- c("bincapt", info.types, "mrds"[is.mrds])

            if (nrow(tem) > 0){
                trap <- tem$trap
                for (i in 1:n){
                    u.id <- unique(id)[i]
                    trig <- trap[id == u.id]
                    out[["bincapt"]][i, trig] <- 1

                    for (j in info.types){
                        for (k in trig){
                            out[[j]][i, k] <- tem[id == u.id & trap == k, j][1]
                        }
                    }

                    if (is.mrds){
                        out[["mrds"]] <- mrds.locs[[s]]
                    }

                }

            }
            out.list[[s]] = out
        }

        if(n.sessions == 1) {
            return(out.list[[1]])
        } else {
            return(out.list)
        }
    } else {
        #if is.animal_ID, output the modified 'captures' directly

        captures$bincapt = 1

        #create a data frame contains all combinations of "session-animal_ID-ID-trap", as captures does not
        #contain the row for the traps without detection for any call
        #and if a session has no detection, then it should appear in the "capture" with n.traps[s] rows with NA "ID"

        tem_session = vector('list', n.sessions)
        for(s in 1:n.sessions){
            if(n_animals[s] > 0){
                tem0 = subset(captures, captures$session == s)
                u.animal = unique(tem0$animal_ID)
                tem_animal = vector('list', n_animals[s])
                for(i in 1:n_animals[s]){
                    animal_label = u.animal[i]
                    tem1 = subset(tem0, tem0$animal_ID == animal_label)
                    u.id = unique(tem1$ID)
                    tem_animal[[i]] = data.frame(session = s, animal_ID = animal_label,
                                                 ID = rep(u.id, each = n.traps[s]),
                                                 trap = rep(1:n.traps[s], length(u.id)))
                }
                tem_session[[s]] = do.call('rbind', tem_animal)
            } else {
                tem_session[[s]] = data.frame(session = s, animal_ID = NA, ID = NA, trap = 1:n.traps[s])
            }
        }

        #this data frame will be naturally well sorted, do not need to sort it again
        head_data = do.call('rbind', tem_session)

        captures = merge(head_data, captures, by = c('session', 'animal_ID', 'ID','trap'), all.x = TRUE)
        captures = captures[,which(colnames(captures) %in% c('session', 'animal_ID', 'ID', 'trap', 'bincapt',
                                                             'bearing', 'toa', 'dist', 'ss', 'mrds_x', 'mrds_y'))]
        for(i in colnames(captures)) captures[,i] = ifelse(is.na(captures[,i]) & !is.na(captures[,'ID']), 0, captures[,i])
        captures = sort.data(captures, 'data.full')
        return(captures)
    }

}



#' A helper function to obtain the distances of the nearest points from a data frame
#'
#' @param from a matrix or a data frame with columns "x" and "y" contains the coordinates of the start points.
#' @param to a matrix or a data frame with columns "x" and "y" contains the coordinates of the end points.
#' @param col_name a character, contains the new column name for the distances in the output, default is
#'                 dist_nearest.
#'
#' @return a data frame as the same as "from", but with an extra column with assigned column name. For each row,
#'         it contains the distance to the nearest point in the "to" data set.
#' @export
#'
#' @examples
dist_nearest = function(from, to, col_name = 'dist_nearest'){
  stopifnot(all(c('x', 'y') %in% colnames(from)))
  stopifnot(all(c('x', 'y') %in% colnames(to)))
  output = as.data.frame(from)
  from = from[, c('x', 'y')]
  to = to[, c('x', 'y')]

  nearest.point = RANN::nn2(to, from, k = 1)

  output[[col_name]] = as.vector(nearest.point$nn.dists)
  return(output)
}





#' Convert traps object
#'
#' Converts an \code{acre} traps matrix to a \code{secr} traps
#' object.
#'
#' The returned object is suitable for use as the \code{traps}
#' argument of the function \link{make.capthist}.
#'
#' @param ss Logical, set to \code{TRUE} if a signal strength
#'     detection function is to be used.
#' @param traps a matrix or a data frame, contains one session's detectors' coordinates
#'
#' @return An object of class \code{traps} comprising a data frame of
#'     x- and y-coordinates, the detector type ('single', 'multi',
#'     'proximity', 'count', 'polygon' etc.), and possibly other
#'     attributes.
#'
#' @export
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

#' Convert mask object
#'
#' Converts an \code{acre} mask matrix to a \code{secr} mask
#' object.
#'
#' The returned object is suitable for use as the \code{mask}
#' argument of the function \link{secr.fit}.
#'
#' @inheritParams fit.acre
#'
#' @return An object of class \code{mask}.
#'
#' @export
convert.mask <- function(mask){
    if (is.list(mask)){
        stop("The convert.mask() function will only convert single-session mask objects.")
    }
    secr::read.mask(data = as.data.frame(mask))
}

#' Capture history conversion.
#'
#' These functions convert a capture history object between the
#' structures required for the \code{acre} and \code{secr}
#' packages.
#'
#' @param capt A \code{secr} capture history object for
#'     \code{convert.capt.to.acre}, or an \code{acre} capture
#'     history object for \code{convert.capt.to.secr}.
#' @param capthist Logical, if \code{TRUE}, a \code{capthist} object
#'     is returned. Otherwise a data frame is returned, which is
#'     suitable for the \code{captures} argument to the
#'     \link{make.capthist} function.
#' @param cutoff The signal strength threshold for detection, if
#'     required.
#' @inheritParams fit.acre
#'
#' @return A capture history object appropriate for analysis using
#'     either the \code{acre} or the \code{secr} package.
#'
#' @name convert.capt
NULL

#' @rdname convert.capt
#' @export
convert.capt.to.acre <- function(capt){
    bincapt <- capt[, 1, ]
    nr <- nrow(bincapt)
    nc <- ncol(bincapt)
    out <- list(bincapt = bincapt)
    if (!is.null(attr(capt, "signalframe"))){
        ss.capt <- matrix(attr(capt, "signalframe")[, 1],
                          nrow = nr, ncol = nc)
        out$ss <- ss.capt
    } else {

    }
    out
}

## Aliasing old convert.capt.to.admbsecr() function name.
#' @rdname convert.capt
#' @export
convert.capt.to.admbsecr <- convert.capt.to.acre

#' @rdname convert.capt
#' @export
convert.capt.to.secr <- function(capt, traps, capthist = TRUE, cutoff = NULL){
    if (!any(names(capt) == "bincapt")){
        if (length(capt) > 1){
            stop("The convert.capt.to.secr() function will only convert single-session capture history objects.")
        }
        capt <- capt[[1]]
    }
    n <- nrow(capt$bincapt)
    n.dets <- sum(capt$bincapt)
    session <- rep(1, n.dets)
    n.indiv.dets <- apply(capt$bincapt, 1, sum)
    ID <- rep(1:n, times = n.indiv.dets)
    occasion <- rep(1, n.dets)
    trap <- c(apply(capt$bincapt, 1, function(x) which(x == 1)),
              recursive = TRUE)
    names(trap) <- NULL
    out <- data.frame(session = session, ID = ID, occasion = occasion,
                      trap = trap)
    fit.ss <- !is.null(capt$ss)
    if (fit.ss){
        out <- data.frame(out, ss = t(capt$ss)[t(capt$bincapt) == 1])
    }
    for (i in names(capt)[!(names(capt) %in% c("bincapt", "ss"))]){
        out[, i] <- t(capt[[i]])[t(capt$bincapt) == 1]
    }
    if (capthist){
        traps <- convert.traps(traps, ss = fit.ss)
        out <- secr::make.capthist(out, traps, fmt = "trapID", noccasions = 1,
                             cutval = cutoff)
    }
    out
}

#' Create a capture history object from a PAMGuard output file
#'
#' Converts a PAMGuard output file to a capture history object
#' suitable for use with the \link{fit.acre} function. This uses
#' \link{make.acoustic.captures} to allocate call identities to
#' detections.
#'
#' @param dets Detection output dataframe from PAMGuard.
#' @param mics A matrix containing the coordinates of microphone
#'     locations.
#' @param time.range A vector of length two, providing a range of
#'     times for which a subset should be taken to create the capture
#'     history.
#' @param sound.speed The speed of sound in metres per second.
#' @param new.allocation Logical, if \code{TRUE}, an improved
#'     call-allocation method is used. The old version is retained so
#'     that older analyses can be replicated.
#' @export
convert.pamguard <- function(dets, mics, time.range = NULL,
                             sound.speed = 330, new.allocation = TRUE){
    mics <- as.matrix(mics)
    toa.info <- dets$startSeconds
    mic.id <- log2(dets$channelMap) + 1
    ss.info <- dets$amplitude
    n <- nrow(dets)
    clicks <- data.frame(session = rep(1, n), ID = 1:n,
                         occasion = rep(1, n), trap = mic.id,
                         ss = ss.info, toa = toa.info)
    if (!is.null(time.range)){
        keep <- toa.info > time.range[1] & toa.info < time.range[2]
        if (!any(keep)){
            stop("No calls were detected within the specified time.range.")
        }
        clicks <- clicks[keep, ]
    }
    ord <- order(clicks$toa)
    clicks <- clicks[ord, ]
    ## Old and new way to allocate IDs below.
    if (new.allocation){
        captures <- clicks
        captures[, 2] <- allocate.calls(mics, clicks, sound.speed)
    } else {
        captures <- make.acoustic.captures(mics, clicks, sound.speed)
    }
    create.capt(captures, traps = mics)
}


#' Create a capture history object from a Raven output file
#'
#' Converts a Raven output file to a capture history object
#' suitable for use with the \link{fit.acre} function. This uses
#' \link{make.acoustic.captures} to allocate call identities to
#' detections.
#'
#' @param dets Detection output dataframe from Raven.
#' @inheritParams convert.pamguard
#' @export
convert.raven <- function(dets, mics, time.range = NULL, sound.speed = 330,
                          new.allocation = TRUE){
    mics <- as.matrix(mics)
    toa.info <- dets[, 4]
    mic.id <- log2(dets[, 3]) + 1
    ss.info <- dets[, 8]
    n <- nrow(dets)
    clicks <- data.frame(session = rep(1, n), ID = 1:n,
                         occasion = rep(1, n), trap = mic.id,
                         ss = ss.info, toa = toa.info)
    if (!is.null(time.range)){
        keep <- toa.info > time.range[1] & toa.info < time.range[2]
        if (!any(keep)){
            stop("No calls were detected within the specified time.range.")
        }
        clicks <- clicks[keep, ]
    }
    ord <- order(clicks$toa)
    clicks <- clicks[ord, ]
    ## Old and new way to allocate IDs below.
    if (new.allocation){
        captures <- clicks
        captures[, 2] <- allocate.calls(mics, clicks, sound.speed)
    } else {
        captures <- make.acoustic.captures(mics, clicks, sound.speed)
    }
    create.capt(captures, traps = mics)
}


#' Assigning ID numbers to sounds
#'
#' Identifies recaptures and assigns ID numbers to sounds recorded for
#' an SECR model.
#'
#' Detected sounds are assumed to come from the same animal if times
#' of arrival at different microphones are closer together than the
#' time it would take for sound to travel between these microphones.
#'
#' @param mics a matrix containing the coordinates of trap locations.
#' @param dets a data frame containing (at least): (i) \code{$toa},
#'     the precise time of arrival of the received sound, and (ii)
#'     \code{$trap} the trap at which the sound was recorded.
#' @param sound.speed the speed of sound in metres per second.
#' @return A data frame. Specifically, the \code{dets} dataframe, now
#'     with a new variable, \code{ID}.
#' @author David Borchers
#'
#' @export
make.acoustic.captures <- function(mics, dets, sound.speed){
    mics <- as.matrix(mics)
    dists <- distances(mics, mics)
    dt <- dists/sound.speed
    K <- dim(mics)[1]
    captures <- dets
    ct <- rep(-Inf, K)
    ID <- 1
    ct[dets$trap[1]] <- dets$toa[1]
    new <- FALSE
    ndets <- length(dets$toa)
    for (i in 2:ndets){
        if (ct[dets$trap[i]] > -Inf){
            nd <- length(which(ct > -Inf))
            captures$ID[(i - nd):(i - 1)] <- ID
            ct <- rep(-Inf, K)
            ct[dets$trap[i]] <- dets$toa[i]
            ID <- ID + 1
            if(i == ndets) captures$ID[i] <- ID
        }
        else {
            ct[dets$trap[i]] <- dets$toa[i]
            ctset <- which(ct > -Inf)
            dts <- dt[ctset, dets$trap[i]]
            cts <- -(ct[ctset] - dets$toa[i])
            if (any((cts - dts) > 0)) new <- TRUE
            if (new) {
                nd <- length(which(ct > -Inf)) - 1
                captures$ID[(i - nd):(i - 1)] <- ID
                ct <- rep(-Inf, K)
                ct[dets$trap[i]] <- dets$toa[i]
                ID <- ID + 1
                new <- FALSE
                if (i == ndets) captures$ID[i] <- ID
            } else if(i == ndets){
                nd <- length(which(ct > -Inf))
                captures$ID[(i - nd + 1):i] <- ID
            }
        }
    }
    captures
}



allocate.calls <- function(mics, dets, sound.speed){
    mics <- as.matrix(mics)
    trap.dists <- distances(mics, mics)
    n.dets <- nrow(dets)
    ## Allocating pairwise plausibility of common cue sources.
    dist.mat <- detection_dists(trap.dists, dets$trap)
    timediff.mat <- detection_timediffs(dets$toa, dets$trap)
    maxtime.mat <- dist.mat/sound.speed
    match.mat <- timediff.mat <= maxtime.mat
    ## Finding blocks of multiple cues with possible common sources.
    incomplete.blocks <- find_incomplete_blocks(match.mat)
    n.blocks <- max(incomplete.blocks)
    complete.block <- logical(n.blocks)
    final.mat <- matrix(FALSE, nrow = n.dets, ncol = n.dets)
    ## Allocating possible common cues to sources.
    reqss.mat <- dist.mat/timediff.mat
    for (i in 1:max(incomplete.blocks)){
        ## Grabbing a block.
        block <- match.mat[incomplete.blocks == i, incomplete.blocks == i]
        reqss <- reqss.mat[incomplete.blocks == i, incomplete.blocks == i]
        ## Working out if there is any possible ambiguity.
        is.complete <- all(block)
        ## If ambiguity, resolve it.
        if (!is.complete){
            block <- blockify(block, reqss)
        }
        final.mat[incomplete.blocks == i, incomplete.blocks == i] <- block
    }
    find_incomplete_blocks(final.mat)
}


convert_one_mask = function(x, trap){
  stopifnot(any(is(x, 'data.frame'), is(x, 'matrix')))

  if(is(x, 'data.frame')){
    x = as.matrix(x)
  }

  stopifnot(ncol(x) == 2)
  stopifnot(is.numeric(x))
  colnames(x) = c('x', 'y')

  if(is.null(attr(x, "buffer"))){
    stopifnot(any(is(trap, 'data.frame'), is(trap, 'matrix')))
    trap = as.matrix(trap)
    stopifnot(is.numeric(trap))
    colnames(trap) = c('x', 'y')

    d = distances(trap, x)
    d_closet_trap = apply(d, 2, min)
    buffer = max(d_closet_trap)
    attr(x, "buffer") = buffer
  }

  if(is.null(attr(x, "area"))){
    sp = as.matrix(dist(x))
    spacing = apply(sp, 1, function(x) min(x[x > 0]))
    spacing = median(spacing, na.rm = T)
    area = spacing^2/10000
    attr(x, "area") = area
  }

  return(x)

}

