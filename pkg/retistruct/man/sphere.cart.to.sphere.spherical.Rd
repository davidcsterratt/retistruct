\name{sphere.cart.to.sphere.spherical}
\alias{sphere.cart.to.sphere.spherical}
\title{Convert from Cartesian to spherical coordinates}
\usage{
  sphere.cart.to.sphere.spherical(P, R = 1)
}
\arguments{
  \item{P}{locations of points on sphere as N-by-3 matrix
  with labelled columns "X", "Y" and "Z"}

  \item{R}{radius of sphere}
}
\value{
  N-by-2 Matrix wtih columns ("phi" and "lambda") of
  locations of points in spherical coordinates
}
\description{
  Convert locations on the surface of a sphere in cartesian
  (X, Y, Z) coordinates to spherical (phi, lambda)
  coordinates.
}
\details{
  It is assumed that all points are lying on the surface of
  a sphere of radius R.
}
\author{
  David Sterratt
}

