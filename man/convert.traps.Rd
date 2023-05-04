% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.r
\name{convert.traps}
\alias{convert.traps}
\title{Convert traps object}
\usage{
convert.traps(traps, ss = FALSE)
}
\arguments{
\item{ss}{Logical, set to \code{TRUE} if a signal strength
detection function is to be used.}
}
\value{
An object of class \code{traps} comprising a data frame of
    x- and y-coordinates, the detector type ('single', 'multi',
    'proximity', 'count', 'polygon' etc.), and possibly other
    attributes.
}
\description{
Converts an \code{ascr} traps matrix to a \code{secr} traps
object.
}
\details{
The returned object is suitable for use as the \code{traps}
argument of the function \link{make.capthist}.
}
\examples{
traps <- convert.traps(traps = example.data$traps)

}