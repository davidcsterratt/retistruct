\name{spherical.to.polar.area}
\alias{spherical.to.polar.area}
\title{Convert lattitude on sphere to radial variable in
area-preserving projection}
\usage{
  spherical.to.polar.area(phi, R = 1)
}
\arguments{
  \item{phi}{Lattitude}

  \item{R}{Radius}
}
\value{
  Coordinate \code{rho} that has the dimensions of length
}
\description{
  Project spherical coordinate system \eqn{(\phi, \lambda)}
  to a polar coordinate system \eqn{(\rho, \lambda)} such
  that the area of each small region is preserved.
}
\details{
  This requires \deqn{R^2\delta\phi\cos\phi\delta\lambda =
  \rho\delta\rho\delta\lambda}.  Hence
  \deqn{R^2\int^{\phi}_{-\pi/2} \cos\phi' d\phi' =
  \int_0^{\rho} \rho' d\rho'}.  Solving gives
  \eqn{\rho^2/2=R^2(\sin\phi+1)} and hence
  \deqn{\rho=R\sqrt{2(\sin\phi+1)}}.

  As a check, consider that total area needs to be
  preserved.  If \eqn{\rho_0} is maximum value of new
  variable then \eqn{A=2\pi
  R^2(\sin(\phi_0)+1)=\pi\rho_0^2}. So
  \eqn{\rho_0=R\sqrt{2(\sin\phi_0+1)}}, which agrees with
  the formula above.
}
\author{
  David Sterratt
}

