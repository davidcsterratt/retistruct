\name{retistruct.export.matlab}
\alias{retistruct.export.matlab}
\title{Save reconstruction data in MATLAB format}
\usage{
retistruct.export.matlab(r)
}
\arguments{
  \item{r}{\code{\link{ReconstructedDataset}} object}
}
\description{
Save as a MATLAB object certain fields of an object
\code{r} that inherits \code{\link{ReconstructedDataset}}
and \code{\link{ReconstructedOutline}} to a file called
\code{r.mat} in the directory \code{r$dataset}.
}
\author{
David Sterratt
}

