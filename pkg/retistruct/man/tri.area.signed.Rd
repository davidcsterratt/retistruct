\name{tri.area.signed}
\alias{tri.area.signed}
\title{"Signed area" of triangles on a plane}
\usage{
  tri.area.signed(P, Pt)
}
\arguments{
  \item{P}{2-column matrix of vertices of triangles}

  \item{Pt}{3-column matrix of indicies of rows of \code{P}
  giving triangulation}
}
\value{
  Vectors of signed areas of triangles. Positive sign
  indicates points are anticlockwise direction; negative
  indicates clockwise.
}
\description{
  "Signed area" of triangles on a plane
}
\author{
  David Sterratt
}
