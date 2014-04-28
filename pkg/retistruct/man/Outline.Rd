\name{Outline}
\alias{Outline}
\title{Outline constructor}
\usage{
Outline(P, scale = NA, im = NULL)
}
\arguments{
  \item{P}{The points of the outline. The last point is not
  repeated.}

  \item{scale}{The length of one unit of \code{P} in
  micrometres. When images are present, this is the length
  of the side of a pixel in the image.}

  \item{im}{An image as a \code{raster} object}
}
\value{
An \code{Outline} object containing the following:
\item{\code{P}}{A N-by-2 matrix of points of the
\code{Outline} arranged in anticlockwise order}
\item{\code{gf}}{For each row of \code{P}, the index of
\code{P} that is next in the outline travelling
anticlockwise (forwards)} \item{\code{gb}}{For each row of
\code{P}, the index of \code{P} that is next in the outline
travelling clockwise (backwards)} \item{\code{im}}{The
image as a \code{raster} object} \item{\code{scale}}{The
length of one unit of \code{P} in micrometres}
}
\description{
Construct an outline object. This sanitises the input
points \code{P}, as described below.
}
\author{
David Sterratt
}

