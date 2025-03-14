% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spheristruct.R
\name{Rcart}
\alias{Rcart}
\title{Restore points to spherical manifold}
\usage{
Rcart(P, R, Rset, i0, phi0, lambda0)
}
\arguments{
\item{P}{Point positions as N-by-3 matrix}

\item{R}{Radius of sphere}

\item{Rset}{Indices of points on rim}

\item{i0}{Index of fixed point}

\item{phi0}{FullCut-off of curtailed sphere in radians}

\item{lambda0}{Longitude of fixed point on rim}
}
\value{
Points projected back onto sphere
}
\description{
Restore points to spherical manifold after an update of the
Lagrange integration rule
}
\author{
David Sterratt
}
