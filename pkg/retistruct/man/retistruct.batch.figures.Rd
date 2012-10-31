\name{retistruct.batch.figures}
\alias{retistruct.batch.figures}
\title{Plot figures for a batch of reconstructions}
\usage{
  retistruct.batch.figures(tldir = ".", outputdir = tldir,
    ...)
}
\arguments{
  \item{tldir}{The top level directory of the tree through
  which to recurse.}

  \item{outputdir}{Directory in which to dump a log file
  and images}

  \item{...}{Parameters passed to plotting functions}
}
\description{
  Recurse through a directory tree, determining whether the
  directory contains valid derived data and plotting graphs
  if it does.
}
\author{
  David Sterratt
}

