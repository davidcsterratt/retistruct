\name{nameLandmark.dataset}
\alias{nameLandmark.dataset}
\title{Name a landmark in a Dataset}
\usage{
  \method{nameLandmark}{dataset} (d, i, name)
}
\arguments{
  \item{d}{\code{dataset} object}

  \item{i}{index of landmark to name}

  \item{name}{name to give landmark}
}
\value{
  New \code{dataset} object in which landmark is named
}
\description{
  Name a landmark in a \code{dataset}. The name of element
  \code{i} of \code{Ss} is set to \code{name}, the name of
  any element that bore the name is set to "" and all other
  elements are unaltered.
}
\author{
  David Sterratt
}

