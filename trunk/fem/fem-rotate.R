## fem-rotate.R - FEM method with rotation 

library(geometry)                       # For trimesh
library(rgl)                       # For 3D plot

## Construct B matrix for a triangle
construct.B <- function(x, y) {
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
construct.D <- function(E=1, nu=0) {
  D <- E * (1-nu^2) * rbind(c(1 , nu, 0       ),
                            c(nu,  1, 0       ),
                            c(0 ,  0, (1-nu)/2))
  return(D)
}

## Construct the stiffness matrix K for a set of triangles defined by
## as set of points P, a triangluation T of those points and the
## elasticity matrix D
##
construct.K <- function(P, T, D) {
  N <- ncol(P)
  M <- nrow(T)
  K <- matrix(0, 2*N, 2*N)
  ## Elasticity matrix for plane stress in an isotropic material
  for (i in 1:M) {
    ## Start of with one triangle
    x <- P[1, T[i,]]
    y <- P[2, T[i,]]

    B <- construct.B(x, y)
    
    inds <- c(2*T[i,1]-1:0,
              2*T[i,2]-1:0,
              2*T[i,3]-1:0)
    K[inds,inds] <- K[inds,inds] + t(B) %*% D %*% B
  }
  return(K)
}

## 2D rotation transformation
rotate <- function(theta) {
  return(rbind(c(cos(theta), -sin(theta)),
               c(sin(theta),  cos(theta))))
}

## 2D shear transformation
shear <- function(k) {
  return(rbind(c(1, k),
               c(0, 1)))
}

## Construct an elasticity matrix
D <- construct.D()

## Construct original set of points
P <- cbind(c(0,   0),
           c(1,   0),
           c(0.5, 1),
           c(1,   1),
           c(0,   1))

## Triangle matrix
T  <- rbind(c(1, 2, 3),
            c(2, 4, 3),
            c(1, 3, 5))
Tt <- rbind(c(1, 2, 3),
            c(2, 4, 3),
            c(1, 3, 4))

## Construct B and K matricies
K <- construct.K(P, T, D)

## Construct new set of points
q1 <- c(0,   0, 0)
q2 <- c(1,   0, 3)
q3 <- c(0.5, 0.2, 2)
q4 <- c(2,   2, 2)
Q <- cbind(q1, q2, q3, q4)

## Compute the force on each vertex in the 3D frame of reference.
##
## P  - points in 2D frame
## T  - triangulation in 2D frame
## K  - stiffness matrix in 2D frame
## Q  - points in 3D frame
## Tt - triangulation in 3D frame
compute.force <- function(P, T, K, Q, Tt) {
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
    ## Hermitian matrix A such that
    ##
    ## M = R %*% A
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
    
    f <- -K[inds, inds] %*% a

    ## Force in frame of new triangle
    g1 <- R %*% f[1:2]
    g2 <- R %*% f[3:4]
    g3 <- R %*% f[5:6]
    g[,Tt[i,]] <- g[,Tt[i,]] + cbind(g1, g2, g3)

    ## Update energy
    E <- E - 0.5*sum(a*f)
  }
  return(list(g=g, E=E))
}

fn <- function(p) {
  Q <- matrix(p, nrow=3)
  o <- compute.force(P, T, K, Q, Tt)
  return(o$E)
}

gr <- function(p) {
  Q <- matrix(p, nrow=3)
  o <- compute.force(P, T, K, Q, Tt)
  return(-o$g)
}

optimise.mapping <- function(Q.init) {
  opt <- optim(as.vector(Q.init), fn, gr, method="BFGS")
  print(opt)
  return(matrix(opt$p, nrow=3))
}

plot.3d <- function(Q, T) {
    rgl.clear()
    triangles3d(matrix(Q[1,t(T)], nrow=3),
                matrix(Q[2,t(T)], nrow=3),
                matrix(Q[3,t(T)], nrow=3),
                color="darkgrey")
}

Q.init <- Q
graphical <- TRUE
for (i in 1:500) {

  o <- compute.force(P, T, K, Q, Tt)
  g <- o$g
  print(o$E)
  
  Q <- Q + 0.01*g

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

