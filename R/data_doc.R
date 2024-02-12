#' Data to showcase a "bearing_dist_hn" demo
#'
#' For the demonstration of the model with bearing and distance as extra information,
#' and half normal as detection function
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "bearing" - extra information, the direction of the location where a call is detected;\cr
#'      "dist" - extra information, the distance to and the location where a call is detected.\cr
#'      }
#'    \item{fix}{a list, contains the coefficient which been fixed instead of estimated by the model}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("bearing_dist_hn")
"bearing_dist_hn"



#' Data to showcase a "bearing_hn" demo
#'
#' For the demonstration of the model with bearing as extra information,
#' and half normal as detection function
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "bearing" - extra information, the direction of the location where a call is detected.\cr
#'      }
#'    \item{fix}{a list, contains the coefficient which been fixed instead of estimated by the model}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("bearing_hn")
"bearing_hn"


#' Data to showcase a "dist_hn" demo
#'
#' For the demonstration of the model with distance as extra information,
#' and half normal as detection function
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "dist" - extra information, the distance to and the location where a call is detected.\cr
#'      }
#'    \item{fix}{a list, contains the coefficient which been fixed instead of estimated by the model}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("dist_hn")
"dist_hn"



#' Data to showcase a "ihd" demo
#'
#' For the demonstration of the model with inhomogeneous density surface,
#' and half normal as detection function
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{par_extend_model}{a list with the name of coefficients to be modeled with additional covariates
#'                            as the elements names, and each element contains a formula for that modeled coefficient.}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{loc.cov}{a data frame contains location related covariates for the modeled coefficients, it must
#'                   contains column 'x' and 'y' as coordinates.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors.\cr
#'      }
#'    \item{fix}{a list, contains the coefficient which been fixed instead of estimated by the model}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ihd")
"ihd"


#' Data to showcase a "ihd_ext" demo
#'
#' For the demonstration of the model with inhomogeneous density surface,
#' and half normal as detection function. One of the coefficients of the detection function, "sigma",
#' is modeled with additional covariates
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{par_extend_model}{a list with the name of coefficients to be modeled with additional covariates
#'                            as the elements names, and each element contains a formula for that modeled coefficient.}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{trap.cov}{a data frame contains detectors related covariates for the modeled coefficients, it must
#'                    contains the indices of detectors, if the detectors are different between any survey sessions,
#'                    the indices of sessions should be provided as well}
#'    \item{loc.cov}{a data frame contains location related covariates for the modeled coefficients, it must
#'                   contains column 'x' and 'y' as coordinates.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors.\cr
#'      }
#'    \item{fix}{a list, contains the coefficient which been fixed instead of estimated by the model}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ihd_ext")
"ihd_ext"






#' Data to showcase a "mul_ses" demo
#'
#' For the demonstration of the model with bearing and distance as extra information,
#' and half normal as detection function. The survey is carried with 2 sessions.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{traps}{a list with a data frame with the coordinates of the acoustic detectors as each of its element}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "bearing" - extra information, the direction of the location where a call is detected;\cr
#'      "dist" - extra information, the distance to and the location where a call is detected.\cr
#'      }
#'    \item{sv}{a list, contains the coefficient which been assigned a start value for modeling}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("mul_ses")
"mul_ses"




#' Data to showcase a "mul_ses_ext" demo
#'
#' For the demonstration of the model with bearing and distance as extra information,
#' and half normal as detection function. Two of the coefficients of the detection function, "g0" and "sigma",
#' are modeled with additional covariates. The survey is carried with 2 sessions.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{par_extend_model}{a list with the name of coefficients to be modeled with additional covariates
#'                            as the elements names, and each element contains a formula for that modeled coefficient.}
#'    \item{traps}{a list with a data frame with the coordinates of the acoustic detectors as each of its element}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{session.cov}{a data frame contains survey sessions related covariates for the modeled coefficients,
#'                       it must contains the indices of sessions}
#'    \item{trap.cov}{a data frame contains detectors related covariates for the modeled coefficients, it must
#'                    contains the indices of detectors, if the detectors are different between any survey sessions,
#'                    the indices of sessions should be provided as well}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "bearing" - extra information, the direction of the location where a call is detected;\cr
#'      "dist" - extra information, the distance to and the location where a call is detected.\cr
#'      }
#'    \item{sv}{a list, contains the coefficient which been assigned a start value for modeling}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("mul_ses_ext")
"mul_ses_ext"


#' Data to showcase a "simple_hhn" demo
#'
#' For the demonstration of the model with hazard half normal detection function.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hhn" - hazard half normal}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("simple_hhn")
"simple_hhn"


#' Data to showcase a "hhn_cue" demo
#'
#' For the demonstration of the model with hazard half normal detection function, and cue rate as additional
#' information to convert the call density to individuals density.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hhn" - hazard half normal}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{survey.length}{a numeric vector or a scalar contains the length of each session.}
#'    \item{cue.rates}{a numeric vector. contains the recorded cue rates in a series of time periods with identical
#'                     length. The length should be equal to the unit length in the "survey.length".}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("hhn_cue")
"hhn_cue"




#' Data to showcase a "simple_hr" demo
#'
#' For the demonstration of the model with hazard rate detection function.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hr" - hazard rate}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors.\cr
#'      }
#'    \item{sv}{a list, contains the coefficient which been assigned a start value for modeling}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("simple_hr")
"simple_hr"


#' Data to showcase a "ss" demo
#'
#' For the demonstration of the model with signal strength detection function, and signal strength
#' as additional information.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "ss" - signal strength}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{ss.opts}{a list contains signal strength model related options. Here it contains "cutoff",
#'                   which indicates the threshold of a distance of 100% detection.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "ss" - the detected signal strength.\cr
#'      }
#'    \item{sv}{a list, contains the coefficient which been assigned a start value for modeling}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ss")
"ss"

#' Data to showcase a "ss_toa" demo
#'
#' For the demonstration of the model with signal strength detection function, and signal strength
#' and time of arriaval as additional information.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "ss" - signal strength}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{ss.opts}{a list contains signal strength model related options. Here it contains "cutoff",
#'                   which indicates the threshold of a distance of 100% detection.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "ss" - the detected signal strength;\cr
#'      "toa" - the time of arrival.\cr
#'      }
#'    \item{sv}{a list, contains the coefficient which been assigned a start value for modeling}
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ss_toa")
"ss_toa"


#' Data to showcase a "ind_bearing_dist" demo
#'
#' For the demonstration of the model with inhomogeneous density surface, bearing and distance as extra information,
#' and half normal as detection function. The coefficients of "alpha", which is used to describe the randomness of distance,
#' is modeled with additional covariates. The survey is carried with 3 sessions, and each records in the capture history
#' contains identity of the detected individual.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hn" - half normal}
#'    \item{par_extend_model}{a list with the name of coefficients to be modeled with additional covariates
#'                            as the elements names, and each element contains a formula for that modeled coefficient.}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{session.cov}{a data frame contains survey sessions related covariates for the modeled coefficients,
#'                       it must contains the indices of sessions}
#'    \item{trap.cov}{a data frame contains detectors related covariates for the modeled coefficients, it must
#'                    contains the indices of detectors, if the detectors are different between any survey sessions,
#'                    the indices of sessions should be provided as well}
#'    \item{survey.length}{a numeric vector or a scalar contains the length of each session.}
#'    \item{control_create_cpat}{a list with arguments for the function "create.capt()" which converts the data input
#'                               to the data to be fed into the model.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "animal_ID" - the indices or identity of the detected individual;\cr
#'      "bearing" - extra information, the direction of the location where a call is detected;\cr
#'      "dist" - extra information, the distance to and the location where a call is detected.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ind_bearing_dist")
"ind_bearing_dist"



#' Data to showcase a "ind_toa_hhn" demo
#'
#' For the demonstration of the model with inhomogeneous density surface, time of arrival as extra information,
#' and hazard half normal as detection function. The survey is carried with 2 sessions, and each records in the
#' capture history contains identity of the detected individual.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "hhn" - hazard half normal}
#'    \item{par_extend_model}{a list with the name of coefficients to be modeled with additional covariates
#'                            as the elements names, and each element contains a formula for that modeled coefficient.}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{loc.cov}{a data frame contains location related covariates for the modeled coefficients, it must
#'                   contains column 'x' and 'y' as coordinates.}
#'    \item{control_create_cpat}{a list with arguments for the function "create.capt()" which converts the data input
#'                               to the data to be fed into the model.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "animal_ID" - the indices or identity of the detected individual;\cr
#'      "toa" - the time of arrival.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ind_toa_hhn")
"ind_toa_hhn"






#' Data to showcase a "ind_ss" demo
#'
#' For the demonstration of the model with inhomogeneous density surface, signal strength as extra information
#' and detection function. One of the coefficients of detection function, "b0.ss", is modeled by additional covariates.
#' The survey is carried with 3 sessions, and each records in the capture history contains identity of the detected individual.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "ss" - signal strength}
#'    \item{par_extend_model}{a list with the name of coefficients to be modeled with additional covariates
#'                            as the elements names, and each element contains a formula for that modeled coefficient.}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{trap.cov}{a data frame contains detectors related covariates for the modeled coefficients, it must
#'                    contains the indices of detectors, if the detectors are different between any survey sessions,
#'                    the indices of sessions should be provided as well}
#'    \item{loc.cov}{a data frame contains location related covariates for the modeled coefficients, it must
#'                   contains column 'x' and 'y' as coordinates.}
#'    \item{ss.opts}{a list contains signal strength model related options. Here it contains "cutoff",
#'                   which indicates the threshold of a distance of 100% detection.}
#'    \item{control_create_cpat}{a list with arguments for the function "create.capt()" which converts the data input
#'                               to the data to be fed into the model.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "animal_ID" - the indices or identity of the detected individual;\cr
#'      "ss" - the detected signal strength.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ind_ss")
"ind_ss"




#' Data to showcase a "ind_ss_log" demo
#'
#' For the demonstration of the model with signal strength as extra information and detection function with the link function of log.
#' The survey is carried with 3 sessions, and each records in the capture history contains identity of the detected individual.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "ss" - signal strength}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{ss.opts}{a list contains signal strength model related options. Here it contains "cutoff",
#'                   which indicates the threshold of a distance of 100% detection, and "ss.link", which
#'                   indicates the link function for the signal strength detection function.}
#'    \item{control_create_cpat}{a list with arguments for the function "create.capt()" which converts the data input
#'                               to the data to be fed into the model.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "animal_ID" - the indices or identity of the detected individual;\cr
#'      "ss" - the detected signal strength.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ind_ss_log")
"ind_ss_log"




#' Data to showcase a "ind_ss_sp" demo
#'
#' For the demonstration of the model with signal strength as extra information and detection function with
#' the link function of spherical. The survey is carried with 3 sessions, and each records in the
#' capture history contains identity of the detected individual.
#'
#' @format a list with all necessary input for the model:
#'  \describe{
#'    \item{detfn}{detechtion function: "ss" - signal strength}
#'    \item{traps}{a data frame with the coordinates of the acoustic detectors}
#'    \item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
#'                               to create the detectable area from the coordinates of the detectors}
#'    \item{ss.opts}{a list contains signal strength model related options. Here it contains "cutoff",
#'                   which indicates the threshold of a distance of 100% detection, and "ss.link", which
#'                   indicates the link function for the signal strength detection function.}
#'    \item{control_create_cpat}{a list with arguments for the function "create.capt()" which converts the data input
#'                               to the data to be fed into the model.}
#'    \item{capt}{a data frame of the capture history, contains columns as follows:\cr
#'      \cr
#'      "session" - the indices of the survey sessions;\cr
#'      "ID" - the indices of the calls been detected;\cr
#'      "occasion" - not been used, could be ignored;\cr
#'      "trap" - the indices of the detectors;\cr
#'      "animal_ID" - the indices or identity of the detected individual;\cr
#'      "ss" - the detected signal strength.\cr
#'      }
#'  }
#'
#' @source created from the simulation
#'
#' @examples
#' o_demo = demo_fit("ind_ss_sp")
"ind_ss_sp"
