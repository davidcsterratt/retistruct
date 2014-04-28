\name{list.datasets}
\alias{list.datasets}
\title{List datasets underneath a directory}
\usage{
list.datasets(path = ".", verbose = FALSE)
}
\arguments{
  \item{path}{Directory path to start searching from}

  \item{verbose}{If \code{TRUE} report on progress}
}
\value{
A vector of directories containing datasets
}
\description{
List valid datasets underneath a directory. This reports
all directories that appear to be valid.
}
\author{
David Sterratt
}

