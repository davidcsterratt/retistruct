\name{azimuthal.equidistant}
\alias{azimuthal.equidistant}
\title{Azimuthal equidistant projection}
\usage{
azimuthal.equidistant(r, ...)
}
\arguments{
  \item{r}{2-column Matrix of spherical coordinates of
  points on sphere. Column names are \code{phi} and
  \code{lambda}.}

  \item{...}{Arguments not used by this projection.}
}
\value{
2-column Matrix of Cartesian coordinates of points on polar
projection. Column names should be \code{x} and \code{y}.
}
\description{
Azimuthal equidistant projection
}
\note{
This is a special case with the point centred on the
projection being the South Pole. The Mathworld equations
are for the more general case.
}
\author{
David Sterratt
}
\references{
\url{http://en.wikipedia.org/wiki/Map_projection},
\url{http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html}
}

