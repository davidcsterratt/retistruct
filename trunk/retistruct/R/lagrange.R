##' Solve the equation of motion 
##'
##' @title Solve Lagrange's (?) equation of motion
##' @param x generalised coordinates of particles
##' @param force generalised force function
##' @param restraint function to restore positions to constraint surface
##' @param mu mass of all particles
##' @param gamma friction coefficient
##' @param dt time step
##' @param Tmax Time at which to end simulations
##' @param record If \code{TRUE} record and return output from all
##' time steps; otherwise just return final values
##' @param ... Extra parameters to give to f()
##' @return list containing generalised coordinates x over time and the time
##' @author David Sterratt
solve.lagrange <- function(x, force, restraint, mu=1, gamma=0, dt, Tmax,
                           record=FALSE) {
  ## Initial conditions
  x0 <- x
  xm <- x

  ## Convenience constant
  epsilon <- gamma/mu/2*dt

  ## Space to record data
  if (record) {
    xs <- matrix(NA, Tmax/dt, length(x))
    ts <- rep(NA, Tmax/dt)
  }
  
  for (i in 1:(Tmax/dt)) {
    if (record) {
      xs[i,] <- x
      ts[i] <- dt*(i-1)
    }
    
    xp <- 1/(1 + epsilon) * (force(x)/mu*dt^2 + 2*x - (1 - epsilon)*xm)
    ## Make sure particles lie on constraint surface
    xp <- restraint(xp)
    xm <- x
    x <- xp
  }

  out <- list(x=x)
  if (record) {
    out <- c(out, list(xs=xs, ts=ts))
  }
  return(out)
}

spring <- function(x) {
  L <- 1
  d <- sqrt((x[1]-x[2])^2)
  return((d-L)*c(x[2]-x[1], x[1]-x[2]))
}
                      
