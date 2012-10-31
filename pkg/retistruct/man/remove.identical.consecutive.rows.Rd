\name{remove.identical.consecutive.rows}
\alias{remove.identical.consecutive.rows}
\title{Remove identical consecutive rows from a matrix}
\usage{
  remove.identical.consecutive.rows(P)
}
\arguments{
  \item{P}{Source matrix}
}
\value{
  Matrix with identical consecutive rows removed.
}
\description{
  This is simlar to unique(), but spares rows which are
  duplicated, but at different points in the matrix
}
\author{
  David Sterratt
}

