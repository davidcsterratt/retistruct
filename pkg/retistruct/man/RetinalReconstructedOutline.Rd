\name{RetinalReconstructedOutline}
\alias{RetinalReconstructedOutline}
\title{RetinalReconstructedOutline constructor}
\usage{
  RetinalReconstructedOutline(r)
}
\arguments{
  \item{r}{Object that inherits
  \code{ReconstructedOutline}}
}
\value{
  \code{RetinalReconstructedOutline} object. This does not
  contain any extra fields, but there are extra mthods
  dthat apply to it.
}
\description{
  Create an object that is specific to retinal datasets.
  This contains methods that return datapoint and landmark
  coordinates that have been transformed according to the
  values of \code{DVflip} and \code{side}.
}
\author{
  David Sterratt
}

