\name{new.ray}
\alias{new.ray}
\title{Compute origin and deflection of new ray}
\usage{
  new.ray(r0, S)
}
\arguments{
  \item{r0}{Vector comprising Cartesian coordinates of
  incident ray origin, angle of incident ray with x-axis
  and index of system surface that origin lies on,
  \code{NA} if the origin is not associated with any
  surface.}

  \item{S}{S System defined by data frame containing left
  margin of surfaces (A), radii of surfaces (R), cut-offs
  (C) and refracive indices (n).}
}
\value{
  Vector comprinsing new origin and angle
}
\description{
  Compute origin and deflection of new ray
}
\author{
  David Sterratt
}

