\name{ReconstructedDataset}
\alias{ReconstructedDataset}
\title{Constructor for RecontructedDataset object}
\usage{
ReconstructedDataset(r, report = message)
}
\arguments{
  \item{r}{Object that of clases
  \code{reconstructedOutline} and \code{dataset}.}

  \item{report}{Function used to report progress.}
}
\value{
\code{\link{ReconstructedDataset}} object containing the
input information and the following modified and extra
information: \item{\code{Dsb}}{Datapoints in barycentric
coordinates} \item{\code{Dsc}}{Datapoints on reconstructed
sphere in cartesian coordinates}
\item{\code{Dss}}{Datapoints on reconstructed sphere in
spherical coordinates} \item{\code{Ssb}}{Landmarks in
barycentric coordinates} \item{\code{Ssc}}{Landmarks on
reconstructed sphere in cartesian coordinates}
\item{\code{Sss}}{Landmarks on reconstructed sphere in
spherical coordinates}
}
\description{
This function infers the coordinates of datapoints \code{Ds
}and landmarks \code{Ss} in spherical coordinates.
}
\author{
David Sterratt
}

