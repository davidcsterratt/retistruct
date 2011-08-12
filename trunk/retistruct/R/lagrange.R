##' Solve the equation of motion 
##'
##' @title Solve Lagrange's (?) equation of motion
##' @param x generalised coordinates of particles
##' @param f generalised force function
##' @param mu mass of all particles
##' @param gamma friction coefficient
##' @param dt time step
##' @param Tmax Time at which to end simulations
##' @param ... Extra parameters to give to f()
##' @return list containing generalised coordinates x over time and the time
##' @author David Sterratt
solve.lagrange <- function(x, f, mu=1, gamma=0, dt, Tmax, ...) {
  ## Initial conditions
  x0 <- x
  xm <- x

  epsilon <- gamma/mu/2*dt
  xs <- matrix(NA, Tmax/dt, length(x))
  ts <- rep(NA, Tmax/dt)
  
  for (i in 1:(Tmax/dt)) {
    xs[i,] <- x
    ts[i] <- dt*(i-1)
    xp <- 1/(1 + epsilon) * (f(x, ...)/mu*dt^2 + 2*x - (1 - epsilon)*xm)
    xm <- x
    x <- xp
  }
  return(list(x=xs, t=ts))
}

spring <- function(x) {
  L <- 1
  d <- sqrt((x[1]-x[2])^2)
  return((d-L)*c(x[2]-x[1], x[1]-x[2]))
}
                      
