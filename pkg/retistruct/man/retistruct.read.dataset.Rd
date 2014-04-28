\name{retistruct.read.dataset}
\alias{retistruct.read.dataset}
\title{Read a retinal dataset}
\usage{
retistruct.read.dataset(dataset, ...)
}
\arguments{
  \item{dataset}{Path to directory containing as SYS and
  MAP file}

  \item{...}{Parameters passed to the format-specific
  functions.}
}
\value{
An object that of classes \code{\link{RetinalDataset}} and
\code{\link{RetinalDataset}}. There may be extra fields
too, depending on the format.
}
\description{
Read a retinal dataset in one of three formats; for
information on formats see see
\code{\link{idt.read.dataset}},
\code{\link{csv.read.dataset}} and
\code{\link{ijroi.read.dataset}}. The format is
autodetected from the files in the directory.
}
\author{
David Sterratt
}

