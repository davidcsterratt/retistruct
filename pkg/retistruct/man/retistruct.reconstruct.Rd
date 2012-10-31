\name{retistruct.reconstruct}
\alias{retistruct.reconstruct}
\title{Reconstruct a retina}
\usage{
  retistruct.reconstruct(o, report = retistruct.report,
    plot.3d = FALSE, dev.flat = NA, dev.polar = NA, ...)
}
\arguments{
  \item{o}{\code{\link{AnnotatedOutline}} object}

  \item{report}{Function to report progress}

  \item{plot.3d}{If \code{TRUE} show progress in a 3D plot}

  \item{dev.flat}{The ID of the device to which to plot the
  flat representation}

  \item{dev.polar}{The ID of the device to which to plot
  the polar representation}

  \item{...}{Parameters to be passed to
  \code{\link{ReconstructedOutline}}}
}
\value{
  Object of classes
  \code{\link{RetinalReconstructedOutline}} and
  \code{\link{RetinalReconstructedDataset}} that contains
  all the reconstruction information
}
\description{
  Reconstruct a retina
}
\author{
  David Sterratt
}

