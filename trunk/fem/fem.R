## E <- 1                                  # Modulus
## nu <- 0                 # Poission's ratio - can be in range (0, 0.5)
construct.K <- function(P, T, E=1, nu=0) {
  N <- nrow(P)
  M <- nrow(T)
  K <- matrix(0, 2*N, 2*N)
  ## Elasticity matrix for plane stress in an isotropic material
  D <- E * (1-nu^2) * rbind(c(1 , nu, 0       ),
                            c(nu,  1, 0       ),
                            c(0 ,  0, (1-nu)/2))
  for (i in 1:M) {
    ## Start of with one triangle
    x <- P[T[i,], 1]
    y <- P[T[i,], 2]

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

    inds <- c(2*T[i,1]-1:0,
              2*T[i,2]-1:0,
              2*T[i,3]-1:0)
    K[inds,inds] <- K[inds,inds] + t(B) %*% D %*% B
  }
  return(K)
}

construct.C <- function(P, T) {
  N <- nrow(P)
  M <- nrow(T)
  C <- matrix(0, 2*N, 2*N)
  for (i in 1:M) {
    ## Start of with one triangle
    x <- P[T[i,], 1]
    y <- P[T[i,], 2]

    C[2*T[i,1]-1:0,2*T[i,2]-1:0] <- diag(2)
    C[2*T[i,2]-1:0,2*T[i,3]-1:0] <- diag(2)
    C[2*T[i,3]-1:0,2*T[i,1]-1:0] <- diag(2)
    C[2*T[i,2]-1:0,2*T[i,1]-1:0] <- diag(2)
    C[2*T[i,3]-1:0,2*T[i,2]-1:0] <- diag(2)
    C[2*T[i,1]-1:0,2*T[i,3]-1:0] <- diag(2)
  }
  return(C)
}

construct.C2 <- function(P, Cu, L) {
  N <- nrow(P)
  M <- length(L)
  C <- matrix(0, 2*N, 2*N)
  for (i in 1:M) {
    C[2*Cu[i,1]-1:0,2*Cu[i,2]-1:0] <- diag(2) / L[i]
    C[2*Cu[i,2]-1:0,2*Cu[i,1]-1:0] <- diag(2) / L[i]
  }
  return(C)
}


rotate <- function(theta) {
  return(rbind(c(cos(theta), -sin(theta)),
               c(sin(theta),  cos(theta))))
}

rotate3 <- function(theta) {
  M <- matrix(0, 6, 6)
  M[1:2, 1:2] <- rotate(theta)
  M[3:4, 3:4] <- rotate(theta)
  M[5:6, 5:6] <- rotate(theta)
  return(M)
}

shear <- function(k) {
  return(rbind(c(1, k),
               c(0, 1)))
}

shear3 <- function(k) {
  M <- matrix(0, 6, 6)
  M[1:2, 1:2] <- shear(k)
  M[3:4, 3:4] <- shear(k)
  M[5:6, 5:6] <- shear(k)
  return(M)
}

solve.fem <- function(K, P, i.fix, P.fix, alpha=1000) {
  ## K %*% a = a
  abar <- matrix(0, nrow(P)*2, 1)
  for (n in 1:length(i.fix)) {
    i <- i.fix[n]
    K[2*i-1:0,2*i-1:0] <- K[2*i-1:0,2*i-1:0] + alpha * diag(2)
    abar[2*i-1:0,1] <- t(P.fix[n,]-P[i,])
  }
  a <- solve(K, alpha * abar)
  Q <- P + matrix(a, nrow(P), 2, byrow=TRUE)
  return(Q)
}

solve.C <- function(C, P, i.fix, P.fix) {
  ind = as.vector(rbind(2*i.fix-1, 2*i.fix))
  print(ind)
  A <- C[-ind, -ind]
  print(dim(A))
  B <- C[-ind,  ind]
  print(dim(B))
  P <- matrix(t(P.fix), ncol(B), 1)
  print(P[1:10,1])
  D <- diag(apply(cbind(A, 2*B), 1, sum))
  print(dim(D))
  print(dim(P))

  Q <- 2 * solve(D - A) %*% B %*% P
  Q <- matrix(Q, nrow(Q)/2, 2, byrow=TRUE)
  R <- matrix(0, nrow(Q) + length(i.fix), 2)
  R[i.fix,] <- P.fix
  R[-i.fix,] <- Q
  return(R)
}

## Points
P <- rbind(c(0, 0),
           c(1, 0),
           c(0.5, 1))
T <- matrix(c(1, 2, 3), 1, 3)

## Original position
v <- matrix(t(P), nrow(P)*2, 1)

K <- construct.K(P, T)

## Translation
u <- v + 1
print("Translation")
a <- u - v                              # Displacement
print("Forces")
print(K %*% a)
print("Strain")
##print(B %*% a)

## Expansion
print("Expansion")
u <- 2 * v
a <- u - v
plot( u[c(1,3,5,1)], u[c(2,4,6,2)], type='l', col="red")
lines(v[c(1,3,5,1)], v[c(2,4,6,2)])
print("Forces")
print(K %*% a)
print("Strain")
##print(B %*% a)

## Rotation
print("Rotation")
u <- rotate3(pi/2) %*% v
a <- u - v
plot( u[c(1,3,5,1)], u[c(2,4,6,2)], type='l', col="red", xlim=c(-1, 1), ylim=c(0,2))
lines(v[c(1,3,5,1)], v[c(2,4,6,2)])
print("Forces")
print(K %*% a)
print("Strain")
##print(B %*% a)

## Shear
print("Shear")
u <- shear3(0.5) %*% v
a <- u - v
plot( u[c(1,3,5,1)], u[c(2,4,6,2)], type='l', col="red", xlim=c(-1, 1), ylim=c(0,2))
lines(v[c(1,3,5,1)], v[c(2,4,6,2)])
print("Forces")
print(K %*% a)
print("Strain")
##print(B %*% a)
