\name{fire}
\alias{fire}
\title{The FIRE algorithm}
\usage{
  fire(r, force, restraint, m = 1, dt = 0.1, maxmove = 100,
    dtmax = 1, Nmin = 5, finc = 1.1, fdec = 0.5,
    astart = 0.1, fa = 0.99, a = 0.1, nstep = 100,
    tol = 1e-05, verbose = FALSE)
}
\arguments{
  \item{r}{Initial locations of particles}

  \item{force}{Force function}

  \item{restraint}{Restraint function}

  \item{m}{Masses of points}

  \item{dt}{Initial time step}

  \item{maxmove}{Maximum distance to move in any time step}

  \item{dtmax}{Maxiumum time step}

  \item{Nmin}{Number of steps after which to start
  increasing \code{dt}}

  \item{finc}{Fractional increase in \code{dt} per time
  step}

  \item{fdec}{Fractional decrease in \code{dt} after a
  stop}

  \item{astart}{Starting value of \code{a} after a stop}

  \item{fa}{Fraction of \code{a} to retain after each step}

  \item{a}{Initial value of \code{a}}

  \item{nstep}{Maxiumum number of steps}

  \item{tol}{Tolerance - if RMS force is below this value,
  stop and report convergence}

  \item{verbose}{If \code{TRUE} report progress verbosely}
}
\value{
  List containing \code{x}, the positions of the points,
  \code{conv}, which is 0 if convergence as occured and 1
  otherwise, and \code{frms}, the root mean square of the
  forces on the particles.
}
\description{
  This is an implmentation of the FIRE algorithm for
  structural relaxation put forward by Bitzek et al. (2006)
}
\author{
  David Sterratt
}
\references{
  Bitzek, E., Koskinen, P., G\"{a}hler, F., Moseler, M.,
  and Gumbsch, P. (2006). Structural relaxation made
  simple. Phys. Rev. Lett., 97:170201.
}
