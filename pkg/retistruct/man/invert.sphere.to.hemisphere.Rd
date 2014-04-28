\name{invert.sphere.to.hemisphere}
\alias{invert.sphere.to.hemisphere}
\title{Invert sphere to hemisphere}
\usage{
invert.sphere.to.hemisphere(r, phi0, ...)
}
\arguments{
  \item{r}{Coordinates of points in spherical coordinates
  represented as 2 column matrix with column names
  \code{phi} (lattitude) and \code{lambda} (longitude).}

  \item{phi0}{The lattitude to map onto the equator}

  \item{...}{Other arguments}
}
\value{
Matrix in same format, but with \code{pi} added to lambda
and \code{phi} negated and scaled so that the longitude
\code{phi0} is projected to 0 degrees (the equator)
}
\description{
Invert image of a partial sphere and scale the longitude so
that points at lattitude \code{phi0} is projected onto a
longitude of 0 degrees (the equator).
}
\author{
David Sterratt
}

