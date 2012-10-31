\name{addClass}
\alias{addClass}
\title{Add new class to class vector of object}
\usage{
  addClass(newclass, obj)
}
\arguments{
  \item{newclass}{New class to add}

  \item{obj}{Object to which to add it}
}
\value{
  New class vector
}
\description{
  This adds a new class to the class vector. If the class
  is dataset type it is prepended at the start of the class
  vector. If it is an outline type, it is put after all the
  dataset classes, but before the other outline ones. This
  is needed for the plotting functions to work properly.
}
\author{
  David Sterratt
}

