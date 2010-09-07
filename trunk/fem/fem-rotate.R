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

construct.D <- function(E=1, nu=0) {
  D <- E * (1-nu^2) * rbind(c(1 , nu, 0       ),
                            c(nu,  1, 0       ),
                            c(0 ,  0, (1-nu)/2))
  return(D)
}

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

rotate <- function(theta) {
  return(rbind(c(cos(theta), -sin(theta)),
               c(sin(theta),  cos(theta))))
}

shear <- function(k) {
  return(rbind(c(1, k),
               c(0, 1)))
}

## Construct an elasticity matrix
D <- construct.D()

## Construct original set of points
p1 <- c(0,   0)
p2 <- c(1,   0)
p3 <- c(0.5, 1)
p4 <- c(1,   1)
P <- cbind(p1, p2, p3, p4)

## Triangle matrix
T <- rbind(c(1, 2, 3),
           c(2, 4, 3))

## Construct B and K matricies
K <- construct.K(P, T, D)

## Construct new set of points
q1 <- c(0,   0, 0)
q2 <- c(1,   0, 3)
q3 <- c(0.5, 0.2, 2)
q4 <- c(2,   2, 2)
Q <- cbind(q1, q2, q3, q4)

compute.force <- function(T, P, Q, K) {
  g <- matrix(0, nrow(Q), ncol(Q))
  for (i in 1:nrow(T)) {
    P0 <- P[,T[i,2:3]] - P[,T[i,1]]
    Q0 <- Q[,T[i,2:3]] - Q[,T[i,1]]

    ## Find transformation M linking two
    ## This means find M in
    ## t(P) %*% t(M) = t(Q)
    ## Which is same as
    ## M %*% P = Q
    M <- t(solve(t(P0), t(Q0)))

    ## Find the rotation matrix R and unitary matrix M such that
    ## M = R %*% U
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

    ##     The SVD decomposition of the matrix as computed by LAPACK/LINPACK,
    ##
    ##                                 X = U D V',                            
    ## According to wikipedia
    ##
    ## R = V D V'

    ## Force in frame of oringinal triangle
    inds <- c(2*T[i,1]-1:0,
              2*T[i,2]-1:0,
              2*T[i,3]-1:0)
    
    f <- -K[inds, inds] %*% a

    ## Force in frame of new triangle
    g1 <- R %*% f[1:2]
    g2 <- R %*% f[3:4]
    g3 <- R %*% f[5:6]
    g[,T[i,]] <- g[,T[i,]] + cbind(g1, g2, g3)
  }
  return(g)
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
for (i in 1:50) {

  g <- compute.force(T, P, Q, K)

  Q <- Q + 0.01*g

  ## print(Q)
  
  if (graphical) {
    plot.3d(Q, T)
    ## plot.default(x=NA, y=NA,
    ##              xlim=c(min(c(P0[1,],P0p[1,],0)), max(c(P0[1,],P0p[1,],0))),
    ##              ylim=c(min(c(P0[2,],P0p[2,],0)), max(c(P0[2,],P0p[2,],0))))
    
    ## trimesh(T, rbind(0, t(P0p)), add=TRUE, col="red")
  }
}
trimesh(T, rbind(0, t(P)))
