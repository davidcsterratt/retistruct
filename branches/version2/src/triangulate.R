if (!require("sp"))  install.packages("sp") # For point.in.polygon
if (!require("geometry")) install.packages("geometry")
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

## Function to return "signed area" of triangles on a plane
## given points P and a triangulation Pt. Positive sign
## indicates points are anticlockwise direction; negative indicates
## clockwise
tri.area.signed <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A, 0)
  BC <- cbind(C-B, 0)
  return(0.5 * extprod3d(AB, BC)[,3])
}

## Function to return area of triangles on a plane
## given points P and a triangulation Pt. Positive sign
## indicates points are anticlockwise direction; negative indicates
## clockwise
tri.area <- function(P, Pt) {
  return(abs(tri.area.signed(P, Pt)))
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

  ## Remove points outwith retinal outline
  Q <- Q[point.in.polygon(Q[,1], Q[,2], P[,1], P[,2])==1,]

  ## Remove points that are within circle wtih diameter given
  ## by the minimum distance between a point and its neighbour
  d1 <- sqrt(apply((P[c(2:nrow(P),1),] - P)^2, 1, sum))
  d2 <- d1[c(nrow(P), 1:(nrow(P)-1))]
  d <- apply(cbind(d1, d2), 1, max)

  ## Distances between each internal point and each vertex
  r2 <- (outer(P[,1], Q[,1], "-"))^2 + (outer(P[,2], Q[,2], "-"))^2 
  Q <- Q[apply(r2 > d^2, 2, all),]
  
  i.et <- 1 ## Indicies of edge triangles
  while (any(i.et)) {
    ## Triangulate
    S <- rbind(P, Q)
    St <- delaunayn(S)

    ## Remove triangles outwith the retinal outline
    Sc <- (S[St[,1],] + S[St[,2],] + S[St[,3],])/3 # Centres
    St <- St[point.in.polygon(Sc[,1], Sc[,2], P[,1], P[,2])==1,]

    ## FIXME: If any triangles have all their verticies in the edge,
    ## split them?? This commented code is a previous attempt
                                        # indicies of edge triangles
    ##   i.et <- which(apply(matrix(St %in% 1:nrow(P), ncol=3), 1, all))
    ##   while(length(i.et) > 0) {
    ##     print("Fixing mesh")
    ##     print(i.et)
    ##     i.new <- nrow(S)+1:length(i.et)
    
    ##     S.new <- (S[St[i.et,1],] + S[St[i.et,2],] + S[St[i.et,3],])/3
    ##     St.new <- rbind(cbind(i.new, St[i.et,c(1,2)]),
    ##                     cbind(i.new, St[i.et,c(2,3)]),
    ##                     cbind(i.new, St[i.et,c(1,3)]))
    ##     St <- St[-i.et,]
    ##     St <- rbind(St, St.new)
    ##     S <- rbind(S, S.new)
    ##     i.et <- which(apply(matrix(St %in% 1:nrow(P), ncol=3), 1, all))
    ##   }
    
    ## Swap orientation of triangles which have clockwise orientation
    areas.signed <- tri.area.signed(S, St)
    St[areas.signed<0,c(2,3)] <- St[areas.signed<0,c(3,2)]
    
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
    A <- Cmat[(nrow(P)+1):nrow(S), (nrow(P)+1):nrow(S)]
    ## i.new
    ## Connections from free points to fixed points
    B <- Cmat[(nrow(P)+1):nrow(S), 1:nrow(P)]
    
    D <- diag(apply(cbind(A, 2 * B), 1, sum))

    ## The matrix Q is an n by 2 matrix comprising the row vectors of the
    ## solution points
    Q <- 2 * solve(D - A) %*% B %*% P

    i.et <- which(apply(matrix(St %in% 1:nrow(P), ncol=3), 1, all))

    Q <- rbind(Q, (S[St[i.et,1],] + S[St[i.et,2],] + S[St[i.et,3],])/3)
    print(i.et)
  }
  
  return(list(St=St, P=P, Q=Q, Cu=Cu, C=C))
}

