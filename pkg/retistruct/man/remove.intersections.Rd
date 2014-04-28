\name{remove.intersections}
\alias{remove.intersections}
\title{Remove intersections between adjacent segements in a closed path}
\usage{
remove.intersections(P, d = 50)
}
\arguments{
  \item{P}{The points, as a 2-column matrix}

  \item{d}{Criterion for maximum distance when points are
  inser}
}
\value{
A new closed path without intersections
}
\description{
Suppose segments AB and CD intersect.  Point B is replaced
by the intersection point, defined B'.  Point C is replaced
by a point C' on the line B'D. The maxium distance of B'C'
is given by the parameter d. If the distance l B'D is less
than 2d, the distance B'C' is l/2.
}
\author{
David Sterratt
}

