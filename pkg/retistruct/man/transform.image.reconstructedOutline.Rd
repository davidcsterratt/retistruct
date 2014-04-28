\name{transform.image.reconstructedOutline}
\alias{transform.image.reconstructedOutline}
\title{Transform an image into the reconstructed space}
\usage{
\method{transform}{image.reconstructedOutline}(r)
}
\arguments{
  \item{r}{\code{reconstructedOutline} object}
}
\value{
\code{reconstructedOutline} object with extra elements
\item{\code{ims}}{Coordinates of corners of pixes in
spherical coordinates} \item{\code{immask}}{Mask matrix
with same dimensions as image \code{im}}
}
\description{
Transform an image into the reconstructed space. The four
corner coordinates of each pixel are transformed into
spherical coordinates and a mask matrix with the same
dimensions as \code{im} is created. This has \code{TRUE}
for pixels that should be displayed and \code{FALSE} for
ones that should not.
}
\author{
David Sterratt
}

