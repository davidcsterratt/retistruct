## dipole.R - quasi van-der-waals force exerted by dipole on point

## scalar product of two column matricies
dot <- function(x, y) {
  return(as.vector(apply(x * y, 1, sum)))
}

## Create a mesh of x and y values of all possible combination of
## the elements of a and b
meshgrid <- function(a,b) {
  list(x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+"))
}

## Compute the dipole potential between a dipole with ends at
## p1 and p2 and a point at q
compute.E <- function(p1, p2, q, d=0.1) {
  ## If there are more points than rows or vice versa, then
  ## replicate the matrix with the smaller number of rows
  if (nrow(q) > nrow(p1)) {
    p1 <- matrix(p1, nrow(q), 2, byrow=TRUE)
    p2 <- matrix(p2, nrow(q), 2, byrow=TRUE)
  }
  if (nrow(p1) > nrow(q)) {
    q <- matrix(q, nrow(p1), 2, byrow=TRUE)
  }
  
  ## Small quantity to add to denominator to prevent NaNs
  epsilon <- 1e-9

  ## Normalised vector along dipole
  v <- (p2 - p1)/sqrt(dot(p2 - p1, p2 - p1))
  s1 <- dot(p1 - q, v)
  s2 <- dot(p2 - q, v)
  r1 <- sqrt(dot(p1 - q, p1 - q))
  r0 <- sqrt((r1+epsilon)^2 - s1^2)

  E <- d/r0*(atan(s2/r0)  - atan(s1/r0)) -
            (asinh(s2/r0) - asinh(s1/r0))

  return(E)
}

## Compute the gradient with respect to q of the dipole potential
## between a dipole with ends at p1 and p2 and a point at q
compute.dEdq <- function(p1, p2, q, d=0.1) {
  if (nrow(q) > nrow(p1)) {
    p1 <- matrix(p1, nrow(q), 2, byrow=TRUE)
    p2 <- matrix(p2, nrow(q), 2, byrow=TRUE)
  }
  if (nrow(p1) > nrow(q)) {
    q <- matrix(q, nrow(p1), 2, byrow=TRUE)
  }
  
  ## Normalised vector along dipole
  v <- (p2 - p1)/sqrt(dot(p2 - p1, p2 - p1))
  s1 <- dot(p1 - q, v)
  s2 <- dot(p2 - q, v)
  p0 <- (s2 * p1 - s1 * p2)/(s2 - s1)
  r0 <- sqrt(dot(p0 - q, p0 - q))
  r1 <- sqrt(dot(p1 - q, p1 - q))
  r2 <- sqrt(dot(p2 - q, p2 - q))

  dEdr0 <- 1/r0 * (-(d/2/r0 * (atan(s2/r0) - atan(s1/r0)) +
                     d/2    * (s2/r2^2     - s1/r1^2))
                   +          (s2/r2       - s1/r1))
  dEds1 <- 1/r1 * (-d/2/r1 + 1)
  dEds2 <- 1/r2 * ( d/2/r2 - 1)
  dEdq <- dEdr0 * (q - p0)/r0 - (dEds1 + dEds2) * v
  
  return(dEdq)
}

## Create three points and two line segments, and then let them move
## with the dipole force
demo.hinge <- function(d=0.01, tmax=1, dt=0.001) {
  ## Connectivity matrix
  Cu <- matrix(c(1, 2, 2, 3), 2, 2, byrow=TRUE)
  ## Desired lengths
  L <- c(1, 1)
  ## Initial position of points
  qs <- matrix(c(-1, 0.1, 0, -0.1, 1, 0.1), 3, 2, byrow=TRUE)

  ## Connectivity matrix
  Cu <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5), 4, 2, byrow=TRUE)
  ## Desired lengths
  L <- c(1, 1, 1, 1)
  ## Initial position of points
  qs <- matrix(c(-2, 0.4, -1, 0.1, 0, -0.1, 1, 0.1, 2, 0.4), 5, 2, byrow=TRUE)


  ## Useful matrices
  C <- rbind(Cu, Cu[,2:1])
  B <- matrix(0, nrow(qs), nrow(C))
  for (i in 1:nrow(C)) {
    B[C[i,1],i] <- 1
  }
  N <- nrow(qs)
  FD <- matrix(0, N, 2)
  
  for (time in 1:tmax) {

    ## Compute dipole force, FD
    ## Find all the line segments interacting with each point
    for (i in 1:N) {
      ## Exclude line segments in which the point occurs
      Cus <- Cu[!apply(Cu == i, 1, any),,drop=FALSE]
      if (nrow(Cus)>0) {
        p1 <- qs[Cus[,1],,drop=FALSE]
        p2 <- qs[Cus[,2],,drop=FALSE]
        q <- qs[i,,drop=FALSE]
        FDs <- -compute.dEdq(p1, p2, q, d)
        FD[i,] <- apply(FDs, 2, sum)
      }
    }
    ## To compensate for the reaction, ensure that the force
    ## on the entire structure is zero 
    FD <- FD - matrix(apply(FD, 2, mean), nrow(FD), 2, byrow=TRUE)

    ## Compute the elastic force, FE
    qi <- qs[C[,1],]
    qj <- qs[C[,2],]
    Deltaq <- (qj - qi)
    l <- sqrt(apply(Deltaq^2, 1, sum))
    fac <- (l - c(L, L))/(c(L, L))
    FE <- 500 * (B %*% (fac * Deltaq))

    ## Plot points and arrows
    plot(qs, xlim=mean(qs[,1]) + c(-2, 2), ylim=mean(qs[,2]) + c(-2, 2))
    text(qs[,1], qs[,2], 1:N)
    
    segments(qs[Cu[,1],1], qs[Cu[,1],2],
             qs[Cu[,2],1], qs[Cu[,2],2])
    arrows(qs[,1], qs[,2], qs[,1] + 0.5 * FD[,1], qs[,2] + 0.5 * FD[,2], len=0.1)
    arrows(qs[,1], qs[,2], qs[,1] + 0.5 * FE[,1], qs[,2] + 0.5 * FE[,2], len=0.1, col="red")    
    ## Update locations using Euler method
    qs <- qs +  dt * (FD + FE)
  }
  print("Desired lengths:")
  print(L)
  print("Actual lengths:")
  print(l[1:(length(l)/2)])
}

## Plot the potential and the gradient due to a diople 
demo.E <- function(d=0.05) {
  xs <- seq(-0.5, 0.5, by=0.05)
  ys <- seq(-0.5, 0.5, by=0.05)
  q.mesh <- meshgrid(xs, ys)
  q <- cbind(as.vector(q.mesh$x), as.vector(q.mesh$y))
  
  ## plot(NA, NA, xlim=c(-1,1), ylim=c(-1, 1))
  ## points(p1[,1], p1[,2])
  ## points(p2[,1], p2[,2])
  ## points(q[,1],  q[,2], col="red")

  ## Create diople with ends at p1 and p2
  p1 <- matrix(c(-0.3, 0), 1, 2)
  p2 <- matrix(c( 0, 0), 1, 2)
  E <- compute.E(p1, p2, q, d)
  gr <- -compute.dEdq(p1, p2, q, d)

  p1 <- matrix(c( 0, 0), 1, 2)
  p2 <- matrix(c( 0.3, 0), 1, 2)
  E  <- E  + compute.E(p1, p2, q, d)
  gr <- gr - compute.dEdq(p1, p2, q, d)

  Emat <- matrix(E, length(xs), length(ys), byrow=TRUE)
  image(xs, ys, Emat ,zlim=c(min(E), 2), col=rainbow(100))

  arrows(q[,1], q[,2], q[,1] + 0.0005 * gr[,1], q[,2] + 0.0005 * gr[,2], len=0.1)
}
