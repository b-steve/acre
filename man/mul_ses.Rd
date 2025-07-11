% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_doc.R
\docType{data}
\name{mul_ses}
\alias{mul_ses}
\title{Data to showcase a "mul_ses" demo}
\format{
a list with all necessary input for the model:
\describe{
\item{detfn}{detechtion function: "hn" - half normal}
\item{traps}{a list with a data frame with the coordinates of the acoustic detectors as each of its element}
\item{control.mask}{a list with the basic argument "buffer", which is the max detectable distance,
to create the detectable area from the coordinates of the detectors}
\item{capt}{a data frame of the capture history, contains columns as follows:\cr
\cr
"session" - the indices of the survey sessions;\cr
"ID" - the indices of the calls been detected;\cr
"occasion" - not been used, could be ignored;\cr
"trap" - the indices of the detectors;\cr
"bearing" - extra information, the direction of the location where a call is detected;\cr
"dist" - extra information, the distance to and the location where a call is detected.\cr
}
\item{sv}{a list, contains the coefficient which been assigned a start value for modeling}
}
}
\source{
created from the simulation
}
\usage{
mul_ses
}
\description{
For the demonstration of the model with bearing and distance as extra information,
and half normal as detection function. The survey is carried with 2 sessions.
}
\keyword{datasets}
