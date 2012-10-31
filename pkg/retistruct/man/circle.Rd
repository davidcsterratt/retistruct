\name{circle}
\alias{circle}
\title{Return points on the unit circle}
\usage{
  circle(n = 12, L = NULL)
}
\arguments{
  \item{n}{Number of points}

  \item{L}{Intervals between points}
}
\value{
  The cartesian coordinates of the points
}
\description{
  Return points on the unit circle in an anti-clockwise
  direction. If \code{L} is not specified \code{n} points
  are returned. If \code{L} is specified, the same number
  of points are returned as there are elements in \code{L},
  the interval between successive points being proportional
  to \code{L}.
}
\author{
  David Sterratt
}

