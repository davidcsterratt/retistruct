\name{sphere.tri.area}
\alias{sphere.tri.area}
\title{Area of triangles on a sphere}
\source{
Wolfram MathWorld
\url{http://mathworld.wolfram.com/SphericalTriangle.html}
and \url{http://mathworld.wolfram.com/SphericalExcess.html}
}
\usage{
sphere.tri.area(P, Pt)
}
\arguments{
  \item{P}{2-column matrix of vertices of triangles given
  in spherical polar coordinates. Columns need to be
  labelled \code{phi} (lattidute) and \code{lambda}
  (longitude).}

  \item{Pt}{3-column matrix of indicies of rows of \code{P}
  giving triangulation}
}
\value{
Vectors of areas of triangles in units of steradians
}
\description{
This uses l'Hullier's theorem to compute the spherical
excess and hence the area of the spherical triangle.
}
\examples{
## Something that should be an eighth of a sphere, i.e. pi/2
P <- cbind(phi=c(0, 0, pi/2), lambda=c(0, pi/2, pi/2))
Pt <- cbind(1, 2, 3)
## The result of this should be 0.5
print(sphere.tri.area(P, Pt)/pi)

## Now a small triangle
P1 <- cbind(phi=c(0, 0, 0.01), lambda=c(0, 0.01, 0.01))
Pt1 <- cbind(1, 2, 3)
## The result of this should approximately 0.01^2/2
print(sphere.tri.area(P, Pt)/(0.01^2/2))

## Now check that it works for both
P <- rbind(P, P1)
Pt <- rbind(1:3, 4:6)
## Should have two components
print(sphere.tri.area(P, Pt))
}
\author{
David Sterratt
}

