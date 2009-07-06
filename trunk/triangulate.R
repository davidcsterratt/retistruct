require("sp")                           # For point.in.polygon
require("geometry")                     # For delauynayn
source("common.R")                      # For meshgrid

connections2triangulation <- function(C) {
  T <- matrix(NA, 0, 3)
  ## Convert a matrix containing links into a triangulation
  for (i in 1:(dim(C)[1])) {
    P1 <- C[i,1]
    P2 <- C[i,2]
    P1n <- setdiff(c(Cu[Cu[,2]==P1, 1], Cu[Cu[,1]==P1, 2]), P2)
    P2n <- setdiff(c(Cu[Cu[,2]==P2, 1], Cu[Cu[,1]==P2, 2]), P1)
##    print(P1n)
##    print(P2n)
    for (P3 in intersect(P1n, P2n)) {
      T.row <- sort(c(P1, P2, P3))
      if (length(unique(T.row)) != 3) {
        print("Error")
        print(T.row)
      } else {
        T <- rbind(T, T.row)
      }
    }
    if (dim(T)[1]==841) {
      print(paste(P1, P2, P3))
      print(paste(P1n))
      print(paste(P2n))
      print(i)
    }

  }
  return(T)
}

## Create grid of N random points around the outline P
create.grid.random <- function(P, N=1000) {
  Q <- cbind(runif(N, min(P[,1]), max(P[,1])),
             runif(N, min(P[,2]), max(P[,2])))

  return(Q)
}

## Create triangular grid of points with separation L around P
create.grid.triangular <- function(P, L=200) {
  Q.mesh <- meshgrid(seq(min(P[,1])/L - 1,             max(P[,1])/L + 1, by=1),
                     seq(min(P[,2])/(sqrt(3)/2*L) - 1, max(P[,2])/(sqrt(3)/2*L) + 1,
                         by=1))
  Q <- cbind(L *            (as.vector(Q.mesh$x) + 0.5 * (as.vector(Q.mesh$y) %% 2)),
             L * sqrt(3)/2 * as.vector(Q.mesh$y))

  return(Q)
}

## Create square grid of points with separation L around P
create.grid.square <- function(P, L) {
  Q.mesh <- meshgrid(seq(min(P[,1])/L - 1, max(P[,1])/L + 1, by=1),
                     seq(min(P[,2])/L - 1, max(P[,2])/L + 1, by=1))
  Q <- cbind(L * as.vector(Q.mesh$x),
             L * as.vector(Q.mesh$y))
  Q <- Q + runif(length(Q), -0.25 * L, 0.25 * L)

  return(Q)
}

## From an outline described by the points P, and a grid creation function,
## create a new set of exteral triangulation points P, a set of internal
## points Q and a triangulation St
create.mesh <- function(P, create.grid=create.grid.random, ...) {
  ## Create grid of points around outline
  Q <- create.grid(P, ...)

  ## Plot them
  ## points(Q, col="red")
                   
  ## Remove points outwith retinal outline
  Q <- Q[point.in.polygon(Q[,1], Q[,2], P[,1], P[,2])==1,]

  ## Remove points that are within circle wtih diameter given
  ## by the minimum distance between a point and its neighbour

  d1 <- sqrt(apply((P[c(2:nrow(P),1),]           - P)^2, 1, sum))
  d2 <- d1[c(nrow(P), 1:(nrow(P)-1))]
  d <- apply(cbind(d1, d2), 1, max)
  
  ## Distances between each internal point and each vertex
  r2 <- (outer(P[,1], Q[,1], "-"))^2 + (outer(P[,2], Q[,2], "-"))^2 
  Q <- Q[apply(r2 > d^2, 2, all),]

  ## Triangulate
  S <- rbind(P, Q)
  St <- delaunayn(S)

  ## Remove triangles outwith the retninal outline
  ## Centres
  Sc <- (S[St[,1],] + S[St[,2],] + S[St[,3],])/3
  St <- St[point.in.polygon(Sc[,1], Sc[,2], P[,1], P[,2])==1,]

  ## Create the asymmetric connectivity matrix
  Cu <- rbind(St[,1:2], St[,2:3], St[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  ## C is the symmetric connectivity matrix
  C <- rbind(Cu, Cu[,2:1])

  ## Create connection matrix
  Cmat <- matrix(0, nrow(S), nrow(S))
  for (i in 1:nrow(C)) {
    Cmat[C[i,1], C[i,2]] <- 1
  }
  ## Connections from free points to free points
  A <- Cmat[nrow(P)+1:nrow(Q), nrow(P)+1:nrow(Q)]
  ## Connections from free points to fixed points
  B <- Cmat[nrow(P)+1:nrow(Q), 1:nrow(P)]
  
  D <- diag(apply(cbind(A, 2 * B), 1, sum))

  ## The matrix Q is an n by 2 matrix comprising the row vectors of the
  ## solution points
  Q <- 2 * solve(D - A) %*% B %*% P

  return(list(St=St, P=P, Q=Q))
}

