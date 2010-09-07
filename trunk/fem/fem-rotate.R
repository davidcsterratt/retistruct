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

rotate <- function(theta) {
  return(rbind(c(cos(theta), -sin(theta)),
               c(sin(theta),  cos(theta))))
}

shear <- function(k) {
  return(rbind(c(1, k),
               c(0, 1)))
}

D <- construct.D()


## Construct original triangle
p1 <- c(0,   0)
p2 <- c(1,   0)
p3 <- c(0.5, 1)
P <- cbind(p2-p1, p3-p1)

P2 <- cbind(p1, p2, p3)
B <- construct.B(P2[1,], P2[2,])
K <- t(B) %*% D %*% B

## Construct new triangle
q1 <- c(0,   0, 0)
q2 <- c(1,   0, 3)
q3 <- c(0.5, 0.2, 2)
Q <- cbind(q2-q1, q3-q1)

## Find transformation M linking two
## This means find M in
## t(P) %*% t(M) = t(Q)
## Which is same as
## M %*% P = Q
M <- t(solve(t(P), t(Q)))

## Find the rotation matrix R and unitary matrix M such that
## M = R %*% U
s <- svd(M)
U <- s$v %*% diag(s$d) %*% t(s$v)
R <- s$u %*% t(s$v)

## Find P', the transformed points in the frame of reference of the
## original triangle
## R %*% Pp = Q
## Hence we need
Pp <- qr.solve(R, Q)

A <- Pp - P
a <- c(0, 0, as.vector(A))

##     The SVD decomposition of the matrix as computed by LAPACK/LINPACK,
##
##                                 X = U D V',                            
## According to wikipedia
##
## R = V D V'

## Force in frame of oringinal triangle

f <- K %*% a

## Force in frame of new triangle
R %*% f[1:2]
