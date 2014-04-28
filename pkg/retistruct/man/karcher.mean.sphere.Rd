\name{karcher.mean.sphere}
\alias{karcher.mean.sphere}
\title{Karcher mean on the sphere}
\usage{
karcher.mean.sphere(x, na.rm = FALSE, var = FALSE)
}
\arguments{
  \item{x}{Matrix of points on sphere as N-by-2 matrix with
  labelled columns \\code{phi} (lattitude) and
  \code{lambda} (longitude)}

  \item{na.rm}{logical value indicating whether \code{NA}
  values should be stripped before the computation
  proceeds.}

  \item{var}{logical value indicating whether variance
  should be returned too.}
}
\value{
Vector of means with components named \code{phi} and
\code{lambda}. If \code{var} is \code{TRUE}, a list
containing mean and varience in elements \code{mean} and
\code{var}.
}
\description{
The Karcher mean of a set of points on a manifold is
defined as the point whose sum of squared Riemmann
distances to the points is minimal. On a sphere using
sphereical coordinates this distance can be computed using
the formula for central angle.
}
\author{
David Sterratt
}
\references{
Heo, G. and Small, C. G. (2006). Form representations and
means for landmarks: A survey and comparative study.
\emph{Computer Vision and Image Understanding},
102:188-203.
}
\seealso{
\code{\link{central.angle}}
}

