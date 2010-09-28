fem.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
}
## fem.R - FEM method with rotation 

## Construct B matrix for a triangle
fem.B <- function(x, y) {
  ## Area of triangle
  Delta <- 0.5 * det(cbind(1, x, y))

  ## Construct B matrix
  a <- c()
  b <- c()
  c <- c()

  a[1] <- x[2]*y[3] - y[2]*x[3]
  a[2] <- x[3]*y[1] - y[3]*x[1]
  a[3] <- x[1]*y[2] - y[1]*x[2]
  b[1] <- y[2] - y[3]
  b[2] <- y[3] - y[1]
  b[3] <- y[1] - y[2]
  c[1] <- x[3] - x[2]
  c[2] <- x[1] - x[3]
  c[3] <- x[2] - x[1]

  B1 <- 1/(2*Delta) * rbind(c(b[1], 0),
                            c(0,    c[1]),
                            c(c[1], b[1]))
  B2 <- 1/(2*Delta) * rbind(c(b[2], 0),
                            c(0,    c[2]),
                            c(c[2], b[2]))
  B3 <- 1/(2*Delta) * rbind(c(b[3], 0),
                            c(0,    c[3]),
                            c(c[3], b[3]))

  B <- cbind(B1, B2, B3)
}

## 
fem.D <- function(E=1, nu=0) {
  D <- E * (1-nu^2) * rbind(c(1 , nu, 0       ),
                            c(nu,  1, 0       ),
                            c(0 ,  0, (1-nu)/2))
  return(D)
}

## Construct the stiffness matrix K for a set of triangles defined by
## as set of points P, a triangluation T of those points and the
## elasticity matrix D
##
fem.K <- function(P, T, D) {
  N <- ncol(P)
  M <- nrow(T)
  K <- matrix(0, 2*N, 2*N)
  ## Elasticity matrix for plane stress in an isotropic material
  for (i in 1:M) {
    ## Start of with one triangle
    x <- P[1, T[i,]]
    y <- P[2, T[i,]]

    B <- fem.B(x, y)
    
    inds <- c(2*T[i,1]-1:0,
              2*T[i,2]-1:0,
              2*T[i,3]-1:0)
    K[inds,inds] <- K[inds,inds] + t(B) %*% D %*% B
  }
  return(K)
}

## Compute the force on each vertex in the 3D frame of reference.
##
## P  - points in 2D frame
## T  - triangulation in 2D frame
## K  - stiffness matrix in 2D frame
## Q  - points in 3D frame
## Tt - triangulation in 3D frame
fem.gradient <- function(P, T, K, Q, Tt) {
  ## Initialise the force vector
  g <- matrix(0, nrow(Q), ncol(Q))
  ## Initialise energy
  E <- 0
  ## Go through each triangle in turn
  for (i in 1:nrow(T)) {
    ## Find the directions of corresponding sides of the triangle in
    ## the flat frame and the 3D frame
    P0 <- P[,T [i,2:3]] - P[,T [i,1]]
    Q0 <- Q[,Tt[i,2:3]] - Q[,Tt[i,1]]

    ## Find transformation M linking two
    ## This means find M in
    ## t(P) %*% t(M) = t(Q)
    ## Which is same as
    ## M %*% P = Q
    M <- t(solve(t(P0), t(Q0)))

    ## Find the unitary rotation matrix R and a positive semidefinite
    ## Hermitian matrix H such that
    ##
    ## M = R %*% H
    ##
    ## To do this, use the svd function, which, for an argument M,
    ## returns a list comprising u, d, v, such that:
    ##
    ## M = U %*% diag(D) %*% t(V)
    ##
    ## The unitary rotation matrix is then given by:
    ##
    ## R = U %*% t(V)
    ##
    ## and the Hermitian matrix is:
    ##
    ## A = V %*% diag(D) %*% t(V)
    ##
    s <- svd(M)
    U <- s$v %*% diag(s$d) %*% t(s$v)
    R <- s$u %*% t(s$v)
    
    ## Find P', the transformed points in the frame of reference of the
    ## original triangle
    ## R %*% Pp = Q
    ## Hence we need
    P0p <- qr.solve(R, Q0)
    
    A <- P0p - P0
    a <- c(0, 0, as.vector(A))

    ## Force in frame of oringinal triangle
    inds <- c(2*T[i,1]-1:0,
              2*T[i,2]-1:0,
              2*T[i,3]-1:0)
    
    f <- K[inds, inds] %*% a

    ## Force in frame of new triangle
    g1 <- R %*% f[1:2]
    g2 <- R %*% f[3:4]
    g3 <- R %*% f[5:6]
    g[,Tt[i,]] <- g[,Tt[i,]] + cbind(g1, g2, g3)

    ## Update energy
    E <- E + 0.5*sum(a*f)
  }
  return(list(g=g, E=E))
}

fem.fn <- function(p) {
  Q <- matrix(p, nrow=3)
  o <- fem.gradient(P, T, K, Q, Tt)
  return(o$E)
}

fem.gr <- function(p) {
  Q <- matrix(p, nrow=3)
  o <- fem.gradient(P, T, K, Q, Tt)
  return(o$g)
}

fem.optimise <- function(P, T, Tt, Q.init) {

  ## Construct an elasticity matrix
  D <- fem.D()
  
  ## Construct B and K matricies
  K <- fem.K(P, T, D)

  opt <- optim(as.vector(Q.init), fem.fn, fem.gr, method="BFGS")
  print(opt)
  return(matrix(opt$p, nrow=3))
}

fem.plot.3d <- function(Q, T) {
    rgl.clear()
    triangles3d(matrix(Q[1,t(T)], nrow=3),
                matrix(Q[2,t(T)], nrow=3),
                matrix(Q[3,t(T)], nrow=3),
                color="darkgrey")
}

fem.solve.euler <- function(P, T, Tt, Q.init,
                            graphical=TRUE) {
  for (i in 1:500) {

    o <- fem.gradient(P, T, K, Q, Tt)
    g <- o$g
    print(o$E)
    
    Q <- Q - 0.01*g

    ## print(Q)
    
    if (graphical) {
      plot.3d(Q, Tt)
      ## plot.default(x=NA, y=NA,
      ##              xlim=c(min(c(P0[1,],P0p[1,],0)), max(c(P0[1,],P0p[1,],0))),
      ##              ylim=c(min(c(P0[2,],P0p[2,],0)), max(c(P0[2,],P0p[2,],0))))
      
      ## trimesh(T, rbind(0, t(P0p)), add=TRUE, col="red")
    }
  }
  trimesh(T, t(P))
  text(P[1,], P[2,], 1:ncol(P))
}

##  Optimisation function when output points are constrained to be on a sphere
fem.optimise.mapping <- function(r, nu=0.45, method="BFGS",
                             plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  phi     <- r$phi
  lambda  <- r$lambda
  R       <- r$R
  phi0    <- r$phi0
  lambda0 <- r$lambda0
  T       <- r$T
  Tt      <- r$Tt
  Rsett   <- r$Rsett
  i0t     <- r$i0t
  P       <- t(r$P)
  Pt      <- r$Pt
  Nt      <- nrow(Pt)  
  Nphi    <- Nt - length(Rsett)


  A       <- r$A
  Cut     <- r$Cut
  Ct      <- r$Ct
  Lt      <- r$Lt
  
  ## Optimisation and plotting
  ## Free up space in Q
  Q <- matrix(0, 3, Nt)

  ## Construct elasticity and stiffness matricies
  D <- fem.D(nu=nu)
  K <- fem.K(P, T, D)

  ## fem.fn.sphere() - requires the following global variables:
  ##
  ## P       - locations of points on plane
  ## T       - triangulation on plane
  ## Tt      - triangulation on sphere
  ## K       - stiffness matrix on plane
  ## R       - radius of sphere
  ## phi0    - lattitude of rim in radians
  ## lambda0 - longitude of fixed point in radians
  ## Nt      - number of points on sphere
  ## Nphi    - number of points on rim
  ## Rsett   - indicies of points sliding on rim
  ## i0t     - index of fixed point on rim
  ##
  fem.fn.sphere <- function(p) {
    ## Extract fixed and non-fixed values of phi and lambda from
    ## variable parameters p
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- p[Nphi+1:(Nt-1)]

    ## Find cartesian coordinates of verticies on sphere
    Q[1,] <- R*cos(phi)*cos(lambda)
    Q[2,] <- R*cos(phi)*sin(lambda)
    Q[3,] <- R*sin(phi)

    ## Find the energy
    o <- fem.gradient(P, T, K, Q, Tt)
    return(o$E)
  }

  fem.gr.sphere <- function(p) {
    ## Extract fixed and non-fixed values of phi and lambda from
    ## variable parameters p
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- p[Nphi+1:(Nt-1)]

    ## Find cartesian coordinates of verticies on sphere
    Q[1,] <- R*cos(phi)*cos(lambda)
    Q[2,] <- R*cos(phi)*sin(lambda)
    Q[3,] <- R*sin(phi)

    ## Find the gradient, with respect to x, y and z
    o <- fem.gradient(P, T, K, Q, Tt)
    dEdQ <- o$g
    
    ## Use the chain rule to compute gradient with respect to phi and lambda
    dQdphi    <- R*rbind(-sin(phi)*cos(lambda),
                         -sin(phi)*sin(lambda),
                         cos(phi))
    dQdlambda <- R*rbind(-cos(phi)*sin(lambda),
                         cos(phi)*cos(lambda),
                         0)
    dEdphi    <- colSums(dQdphi    * dEdQ)
    dEdlambda <- colSums(dQdlambda * dEdQ)
    
    ## Return only the variable parts of this
    return(c(dEdphi[-Rsett], dEdlambda[-i0t]))
  }

  fem.report <- function() {
    ## Report
    print(paste("Strain energy:", fem.fn.sphere(opt$p)))

    ## Report on flipped triangles
    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    print(paste(nflip, "flipped triangles:"))
    print(which(ft$flipped))
    print("Areas")
    print(ft$areas[ft$flipped])
    print(A[ft$flipped])
    print(length(A))
  }
  
  opt <- list()
  opt$p <- c(phi[-Rsett], lambda[-i0t])
  opt$conv <- 1

  ## Report
  fem.report()
  
  while (opt$conv) {
    ## Optimise
    opt <- optim(opt$p, fn=fem.fn.sphere, gr=fem.gr.sphere,
                 method=method)
    
    ## Decode p vector
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- opt$p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

    ## Report
    fem.report()
    
    ## Plot
    if (plot.3d) {
      plot.sphere.spherical(phi, lambda, R, Tt, Rsett)
      plot.outline.spherical(phi, lambda, R, r$gb, r$ht)
    }

    if (!is.na(dev.grid)) {
      dev.set(dev.grid)
      with(r, plot.outline.flat(P, gb))
      plot.gridlines.flat(r$P, r$T, phi, lambda, Tt, phi0*180/pi)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      plot.polar(phi0 * 180/pi)
      r$phi <- phi
      r$lambda <- lambda
      plot.outline.polar(r)
    }
  }
  ft <- flipped.triangles(phi, lambda, Tt, R)
  return(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped)))
}
