\name{retistruct.read.markup}
\alias{retistruct.read.markup}
\title{Read the markup data}
\usage{
  retistruct.read.markup(a, error = stop)
}
\arguments{
  \item{a}{Dataset object, containing \code{dataset} path}

  \item{error}{Function to run on error, by default
  \code{stop()}}
}
\value{
  o \code{RetinalDataset} object \item{V0}{Indicies in
  \code{P} of apicies of tears} \item{VB}{Indicies in
  \code{P} of backward verticies of tears}
  \item{VF}{Indicies in \code{P} of backward verticies of
  tears} \item{iN}{Index in \code{P} of nasal point, or
  \code{NA} if not marked} \item{iD}{Index in \code{P} of
  dorsal point, or \code{NA} if not marked}
  \item{iOD}{Index in \code{Ss} of optic disc }
  \item{phi0}{Angle of rim in degrees}
  \item{DVflip}{Boolean variable indicating if DV axis has
  been flipped}
}
\description{
  Read the markup data contained in the files
  \file{markup.csv}, \file{P.csv} and \file{T.csv} in the
  directory \file{dataset}, which is specified in the
  reconstruction object \code{r}.
}
\details{
  The tear information is contained in the files
  \file{P.csv} and \file{T.csv}. The first file contains
  the locations of outline points that the tears were
  marked up on. The second file contains the indicies of
  the apicies and backward and forward verticies of each
  tear. It is necessary to have the file of points just in
  case the algorithm that determines \code{P} in
  \code{\link{retistruct.read.dataset}} has changed since
  the markup of the tears.

  The remaining information is contained containted in the
  file \file{markup.csv}.

  If \code{DVflip} is specified, the locations of points
  \code{P} flipped in the \eqn{y}-direction. This operation
  also requires the swapping of \code{gf} and \code{gb} and
  \code{VF} and \code{VB}.
}
\author{
  David Sterratt
}

