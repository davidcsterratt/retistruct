fire <- function(r, force, restraint, dt=0.1, maxmove=1E2, dtmax=10,
                 Nmin=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99, a=0.1,
                 nstep=100) {
  Nsteps <- 0 
  # Initialise velocity
  v <- 0*r                              

  for (i in 1:nstep) {
    f <- force(r)
    vf <- dot(f, v)
    ## if (all(vf > 0)) {
    ## vf <- sum(f*v)
    if (sum(vf) > 0) {
      v[vf<0,] <- 0
      v <- (1 - a)*v + a*f/vecnorm(f)*vecnorm(v)
      if (Nsteps > Nmin) {
        dt <- min(dt*finc, dtmax)
        if (dt==dtmax)
          message("Hitting time limit")
        a <- a*fa
      }
      Nsteps <- Nsteps + 1
    } else {
      message("Stop")
      v[vf<0,] <- 0
      a <- astart
      dt <- dt*fdec
      Nsteps <- 0
    }
    v <- v + dt*f
    dr <- dt*v
    normdr <- sqrt(sum(dr*dr))
    if (normdr > maxmove) {
      message("Hitting maxmove limit")
      dr <- maxmove*dr/normdr
    }
    r <- r + dr
    r <- restraint(r)
  }
  return(list(x=r))
}
