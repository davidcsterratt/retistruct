\name{ray.arc.intersection}
\alias{ray.arc.intersection}
\title{Compute intersection of ray and arc}
\usage{
  ray.arc.intersection(Ptheta, A, R, C)
}
\arguments{
  \item{Ptheta}{Vector comprising cartesian coordinates of
  incident ray origin and angle with x-axis of incident
  ray}

  \item{A}{Left-edge of circle on x-axis}

  \item{R}{Radius of circle}

  \item{C}{Angle at which arc cut off}
}
\value{
  Vector of x, y and lambda
}
\description{
  Compute intersection of ray and arc
}
\author{
  David Sterratt
}

