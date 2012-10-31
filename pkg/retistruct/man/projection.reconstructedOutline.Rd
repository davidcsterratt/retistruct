\name{projection.reconstructedOutline}
\alias{projection.reconstructedOutline}
\title{Projection of a reconstructed outline}
\usage{
  \method{projection}{reconstructedOutline} (r,
    transform = identity.transform,
    axisdir = cbind(phi = 90, lambda = 0),
    projection = azimuthal.equalarea,
    proj.centre = cbind(phi = 0, lambda = 0),
    lambdalim = c(-180, 180), philim = c(-90, 90),
    labels = c(0, 90, 180, 270), grid = TRUE,
    grid.col = "gray", grid.bg = "transparent",
    grid.int.minor = 15, grid.int.major = 45,
    colatitude = TRUE, pole = FALSE, image = TRUE, ...)
}
\arguments{
  \item{r}{\code{ReconstructedOutline} object}

  \item{transform}{Transform function to apply to spherical
  coordinates before rotation}

  \item{axisdir}{Direction of axis (North pole) of sphere
  in external space as matrix with column names \code{phi}
  (elevation) and \code{lambda} (longitude).}

  \item{projection}{Projection in which to display object,
  e.g. \code{\link{azimuthal.equalarea}} or
  \code{\link{sinusoidal}}}

  \item{proj.centre}{Location of centre of projection as
  matrix with column names \code{phi} (elevation) and
  \code{lambda} (longitude).}

  \item{lambdalim}{Limits of longitude (in degrees) to
  display}

  \item{philim}{Limits of latitude (in degrees) to display}

  \item{labels}{Vector of 4 labels to plot at 0, 90, 180
  and 270 degrees}

  \item{grid}{Whether or not to show the grid lines of
  lattitude and longitude}

  \item{grid.col}{Colour of the minor grid lines}

  \item{grid.bg}{Background colour of the grid}

  \item{grid.int.minor}{Interval between minor grid lines
  in degrees}

  \item{grid.int.major}{Interval between major grid lines
  in degrees}

  \item{colatitude}{If \code{TRUE} have radial labels
  plotted with respect to colatitude rather than latitude}

  \item{pole}{If \code{TRUE} indicate the pole with a "*"}

  \item{image}{If \code{TRUE}, show the image}

  \item{...}{Graphical parameters to pass to plotting
  functions}
}
\description{
  Draw a projection of a
  \code{\link{ReconstructedOutline}}. This method sets up
  the grid lines and the angular labels and draws the
  image.
}

