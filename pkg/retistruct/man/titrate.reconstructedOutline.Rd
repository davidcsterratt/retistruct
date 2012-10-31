\name{titrate.reconstructedOutline}
\alias{titrate.reconstructedOutline}
\title{Titrate values of phi0}
\usage{
  titrate.reconstructedOutline(r, alpha = 8, x0 = 0.5,
    byd = 1, len.up = 5, len.down = 20)
}
\arguments{
  \item{r}{\code{\link{ReconstructedOutline}} object}

  \item{alpha}{Area penalty scaling coefficient}

  \item{x0}{Area cutoff coefficient}

  \item{byd}{Increments in degrees}

  \item{len.up}{How many increments to go up from starting
  value of \code{phi0} in \code{r}.}

  \item{len.down}{How many increments to go up from
  starting value of \code{phi0} in \code{r}.}
}
\value{
  dat Output data frame
}
\description{
  Try a range of values of phi0s in the reconstruction,
  recording the energy of the mapping in each case.
}
\author{
  David Sterratt
}

