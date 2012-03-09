## Kernel density clustering 

##' @title Kernel density estimate on sphere using Fisherian density
##' with Cartesian coordinates
##' @param r Locations at which to estimate density in Cartesian
##' coordinates on unit sphere
##' @param mu Locations of data points in Cartesian coordinates on
##' unit sphere
##' @param kappa Concentration parameter
##' @return Vector of density estimates
##' @author David Sterratt
##' @export
kde.fhat.cart <- function(r, mu, kappa) {
  if (is.vector(r)) {
    if (length(r) == 3) {
      r <- matrix(r, ncol=3)
    } else {
      stop("r does not have 3 elements")
    }
  } else {
    if (ncol(r) != 3) {
      stop("r does not have 3 columns")
    }
  }

  n <- nrow(mu)
  if (kappa < 1e-10) {
    fac <- 1
  } else {
    fac <- kappa/sinh(kappa)
  }
  return(fac/(4*pi*n)*rowSums(exp(kappa*(outer(r[,1], mu[,1]) +
                                         outer(r[,2], mu[,2]) +
                                         outer(r[,3], mu[,3])))))
}

##' @title Kernel density estimate on sphere using Fisherian density
##' with polar coordinates
##' @param r Locations at which to estimate density in polar
##' coordinates
##' @param mu Locations of data points in polar coordinates
##' @param kappa Concentration parameter
##' @return Vector of density estimates
##' @author David Sterratt
##' @export
kde.fhat <- function(r, mu, kappa) {
  return(kde.fhat.cart(sphere.spherical.to.sphere.cart(r[,"phi"],  r[,"lambda"]),
                       sphere.spherical.to.sphere.cart(mu[,"phi"], mu[,"lambda"]),
                       kappa))
}


##' @title Estimate of the log likelihood of the points mu given a
##' particular value of the concentration kappa
##' @param mu Locations of data points in Cartesian coordinates on
##' unit sphere
##' @param kappa Concentration parameter
##' @return Log likelihood of data
##' @author David Sterratt
##' @export
kde.L <- function(mu, kappa) {
  ## We get clever here, and define a funtion within a function. This
  ## is the kernel density for data point i if the density is
  ## determined by all the other points. Note that mu[-i,] means all
  ## the rows of mu apart from row i
  log.fhati <- function(i) {
    return(log(kde.fhat.cart(mu[i,], mu[-i,], kappa)))
  }

  ## We now pass this function to sapply, which creates a vector of
  ## the result of log.fhati for every value of i
  return(sum(sapply(1:nrow(mu), log.fhati)))
}

##' @title Find the optimal concentration for a set of data
##' @param mu Data in spherical coordinates
##' @return The optimal concentration
##' @author David Sterratt
##' @export
kde.compute.concentration <- function(mu) {
  mu.cart <- sphere.spherical.to.sphere.cart(mu[,"phi"], mu[,"lambda"]) 
  opt <- optimise(function(kappa) {kde.L(mu.cart, kappa)}, interval=c(0, 500),
                  maximum=TRUE)
  return(opt$maximum)
}





