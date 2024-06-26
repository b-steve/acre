% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{distances}
\alias{distances}
\title{Calculating distances between two sets of points}
\usage{
distances(a, b)
}
\arguments{
\item{a}{matrix containing a set of coordinates.}

\item{b}{matrix containing another set of coordinates.}
}
\value{
A matrix with pairwise distances between the two sets of points.
}
\description{
Calculates pairwise distances between two sets of points.
}
\examples{
dists <- distances(example.data$traps, example.data$mask)

}
