fire <- function(r, force, restraint, dt=0.1, maxmove=0.2, dtmax=1.0,
                 Nmin=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99, a=0.1,
                 nstep=100) {
  Nsteps <- 0 
  # Initialise velocity
  v <- 0*r                              

  for (i in 1:nstep) {
    f <- force(r)
    vf <- sum(f*v)
    if (vf > 0.0) {
      v <- (1.0 - a)*v + a*f/vecnorm(f)*vecnorm(v)
      if (Nsteps > Nmin) {
        dt <- min(dt*finc, dtmax)
        a <- a*fa
      }
      Nsteps <- Nsteps + 1
    } else {
      v <- v*0
      a <- astart
      dt <- dt*fdec
      Nsteps <- 0
    }
    v <- v + dt*f
    dr <- dt*v
    normdr <- sqrt(sum(dr*dr))
    if (normdr > maxmove) {
      dr <- maxmove*dr/normdr
    }
    r <- r + dr
    r <- restraint(r)
  }
  return(list(x=r))
}
