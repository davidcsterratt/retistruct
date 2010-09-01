## Now for the dreaded elastic error function....
E2 <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
              E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  cosp <- cos(phi)
  sinp <- sin(phi)
  cosl <- cos(lambda)
  sinl <- sin(lambda)

  x    <-  R*cosp*cosl
  y    <-  R*cosp*sinl
  z    <-  R*sinp

  E.E <- 0
  for (i in 1:nrow(Cu)) {
    j <- Cu[i,1]
    k <- Cu[i,2]
    l <- sqrt((x[j]-x[k])^2 + 
              (y[j]-y[k])^2 +
              (z[j]-z[k])^2)
    E.E <- E.E + 0.5 * (l - L[i])^2/L[i]
  }

  E.A <- 0
  if (E0.A) {
  }
  return(E.E + E0.A*E.A)
}


## Alternative, sequential version of the gradient
dE2 <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
                E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]
  
  cosp <- cos(phi)
  sinp <- sin(phi)
  cosl <- cos(lambda)
  sinl <- sin(lambda)

  x    <-  R*cosp*cosl
  y    <-  R*cosp*sinl
  z    <-  R*sinp

  dxdp <- -R*sinp*cosl
  dydp <- -R*sinp*sinl
  dzdp <-  R*cosp

  dxdl <- -R*cosp*sinl
  dydl <-  R*cosp*cosl
  dzdl <-  rep(0, N)

  dEdp <- rep(0, N)
  dEdl <- rep(0, N)

  if (verbose) print(paste(x[10], y[10], z[10]))
  
  for (i in 1:nrow(Cu)) {
    j <- Cu[i,1]
    k <- Cu[i,2]
    # print(paste("i", i))
    # print(j)
    # print(k)

    dx <- x[j]-x[k]
    dy <- y[j]-y[k]
    dz <- z[j]-z[k]
    
    l <- sqrt(dx^2 + dy^2 + dz^2)
    
    dlidpj <-  1/l*(dx*dxdp[j] + dy*dydp[j] + dz*dzdp[j])
    dlidpk <- -1/l*(dx*dxdp[k] + dy*dydp[k] + dz*dzdp[k])
    dlidlj <-  1/l*(dx*dxdl[j] + dy*dydl[j] + dz*dzdl[j])
    dlidlk <- -1/l*(dx*dxdl[k] + dy*dydl[k] + dz*dzdl[k])

    #print(dlidpj)
#    print(dlidpk)
#    print(dlidlj)
#    print(dxdl[j])
#    print(dydl[j])
#    print(dzdl[j])
  #  print(dlidlk)
    dll <- (l - L[i])/L[i]
      
    dEdp[j] <- dEdp[j] + dlidpj * dll
    dEdp[k] <- dEdp[k] + dlidpk * dll
    dEdl[j] <- dEdl[j] + dlidlj * dll
    dEdl[k] <- dEdl[k] + dlidlk * dll
  }
#  print(dEdl)
  dEAdp <- rep(0, N)
  dEAdl <- rep(0, N)


  if (E0.A) {
  }
  
  return(c(dEdp[-Rset]  + E0.A * dEAdp[-Rset],
           dEdl[-i0]    + E0.A * dEAdl[-i0]))
  ## return(c(dEdp, dEdl))
}

## Alternative, sequential version of the gradient
dE2c <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
                E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  dE <- .Call("spheristruct_E2",
              phi, lambda, Cu, L, R)

  dEdp <- dE[1:N]
  dEdl <- dE[N+(1:N)]
  return(c(dEdp[-Rset], dEdl[-i0]))
  #return(dE[-c(Rset, i0+Nphi)])
  ## return(dE)
}
