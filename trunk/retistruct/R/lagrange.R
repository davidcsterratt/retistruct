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

solve.lagrange.adaptive <- function(x, force, restraint, mu=1, gamma=0, dt,
                                    Delta0=1e-3, nstep) {
  ## Initial conditions
  x0  <- x
  xm  <- x
  xm2 <- x
  t <- 0

  for (i in 1:nstep) {
    ## Convenience constant
    epsilon <- gamma/mu/2*dt

    ## Compute acceleration
    a <- force(x)/mu
    
    ## Assess accuracy

    ## Take two single steps 
    xp <- 1/(1 +   epsilon)*(a*dt^2     + 2*x  - (1 -   epsilon)*xm )
    xp <- restraint(xp)
    
    x2 <- 1/(1 +   epsilon)*(a*dt^2     + 2*xp - (1 -   epsilon)*x  )
    x2 <- restraint(x2)
    
    ## Take one double step
    x1 <- 1/(1 + 2*epsilon)*(a*(2*dt)^2 + 2*x  - (1 - 2*epsilon)*xm2)
    x1 <- restraint(x1)

    ## Compute error estimate
    Delta <- max(abs(x2 - x1))

    ## Compute optimal dt
    dt0 <- ((Delta0/Delta)^0.25)*dt

    ## If dt0 is at least twice the size of dt, we double dt
    if (dt0 > 2*dt) {
      dt <- 2*dt
      xm <- x
      ## xm2 <- xm2
      x <- x1
      t <- t + dt
      message(paste("Doubling time step to ", dt))
    } else {  
      ## If dt0 is less than dt, we half dt
      if (dt0 < dt) {
        dt <- dt/2
        xm2 <- xm
        xm <- 0.5*(x + xm) - a*dt^2
        message(paste("Halving time step to ", dt))
      } else {
        ## Keep going
        xm2 <- xm
        xm <- x
        x <- xp
        t <- t + dt
      }
    }
  }

  out <- list(x=x)
  return(out)
}


spring <- function(x) {
  L <- 1
  d <- sqrt((x[1]-x[2])^2)
  return((d-L)*c(x[2]-x[1], x[1]-x[2]))
}
                      
