\name{sphere.wedge.to.sphere.cart}
\alias{sphere.wedge.to.sphere.cart}
\title{Convert from 'wedge' to Cartesian coordinates}
\usage{
sphere.wedge.to.sphere.cart(psi, f, phi0, R = 1)
}
\arguments{
  \item{psi}{vector of slice angles of N points}

  \item{f}{vector of fractional distances of N points}

  \item{phi0}{rim angle as colatitude}

  \item{R}{radius of sphere}
}
\value{
An N-by-3 matrix in which each row is the cartesian (X, Y,
Z) coordinates of each point
}
\description{
This in the inverse of
\code{\link{sphere.cart.to.sphere.wedge}}
}
\author{
David Sterratt
}

