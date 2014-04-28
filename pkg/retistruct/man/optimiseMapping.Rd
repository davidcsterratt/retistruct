\name{optimiseMapping}
\alias{optimiseMapping}
\title{Optimise mapping}
\usage{
optimiseMapping(r, alpha = 4, x0 = 0.5, nu = 1, method = "BFGS",
  plot.3d = FALSE, dev.flat = NA, dev.polar = NA, control = list())
}
\arguments{
  \item{r}{reconstructedOutline object}

  \item{alpha}{Area penalty scaling coefficient}

  \item{x0}{Area penalty cut-off coefficient}

  \item{nu}{Power to which to raise area}

  \item{method}{Method to pass to \code{optim}}

  \item{plot.3d}{If \code{TRUE} make a 3D plot in an RGL
  window}

  \item{dev.flat}{Device handle for plotting flatplot
  updates to. If \code{NA} don't make any flat plots}

  \item{dev.polar}{Device handle for plotting polar plot
  updates to. If \code{NA} don't make any polar plots.}

  \item{control}{Control argument to pass to \code{optim}}
}
\value{
reconstructedOutline object
}
\description{
Optimise the mapping from the flat outline to the sphere
}
\author{
David Sterratt
}

