\name{retistruct.batch}
\alias{retistruct.batch}
\title{Batch operation using the multicore package}
\usage{
  retistruct.batch(tldir = ".", outputdir = tldir,
    datasets = NULL, device = "pdf", titrate = FALSE,
    cpu.time.limit = 3600, mc.cores = getOption("cores"))
}
\arguments{
  \item{tldir}{If datasets is not specified, the top level
  of the directory tree through which to recurse in order
  to find datasets.}

  \item{outputdir}{directory in which to dump a log file
  and images}

  \item{datasets}{Vector of dataset directories to
  reconstruct}

  \item{device}{string indicating what type of graphics
  output required. Options are "pdf" and "png".}

  \item{titrate}{Whether to "titrate" the reconstruction
  for different values of \code{phi0}. See
  \code{titrate.reconstructedOutline}.}

  \item{cpu.time.limit}{amount of CPU after which to
  terminate the process}

  \item{mc.cores}{The number of cores to use. Defaults to
  the total number available.}
}
\description{
  This function reconstructs a number of datasets, using
  the R \code{multicore} package to distribute the
  reconstruction of multiple datasets across CPUs. If
  \code{datasets} is not specified the function recurses
  through a directory tree starting at \code{tldir},
  determining whether the directory contains valid raw data
  and markup, and performing the reconstruction if it does.
}
\author{
  David Sterratt
}
