\name{flipped.triangles}
\alias{flipped.triangles}
\title{Determine indicies of triangles that are flipped}
\usage{
  flipped.triangles(phi, lambda, Tt, R)
}
\arguments{
  \item{phi}{Vector of lattitudes of points}

  \item{lambda}{Vector of longitudes of points}

  \item{Tt}{Triangulation of points}

  \item{R}{Radius of sphere}
}
\value{
  List containing: \item{\code{flipped}}{Indicies of in
  rows of \code{Tt} of flipped triangles.}
  \item{\code{cents}}{Vectors of centres.}
  \item{\code{areas}}{Areas of triangles.}
}
\description{
  In the projection of points onto the sphere, some
  triangles maybe flipped, i.e. in the wrong orientation.
  This functions determines which triangles are flipped by
  computing the vector pointing to the centre of each
  triangle and comparing this direction to vector product
  of two sides of the triangle.
}
\author{
  David Sterratt
}
