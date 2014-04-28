\name{retistruct.batch.export.matlab}
\alias{retistruct.batch.export.matlab}
\title{Export data from reconstruction data files to matlab}
\usage{
retistruct.batch.export.matlab(tldir = ".")
}
\arguments{
  \item{tldir}{The top level of the directory tree through
  which to recurse}
}
\description{
Recurse through a directory tree, determining whether the
directory contains valid derived data and converting
r.rData files to files in matlab format named r.mat
}
\author{
David Sterratt
}

