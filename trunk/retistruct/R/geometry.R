## scalar product of two column matricies
dot <- function(x, y) {
  return(rowSums(x * y))
}

## Distance norm
norm <- function(X) {
  if (is.vector(X)) {
    return(sqrt(sum(X^2)))
  } else {
    return(sqrt(rowSums(X^2)))
  }
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

## Check whether line from P1 to Q1 and from P2 to Q2 intersect
find.intersection <- function(P1, Q1, P2, Q2) {
  if ((max(P1[1], Q1[1]) < min(P2[1], Q2[1])) ||
      (min(P1[1], Q1[1]) > max(P2[1], Q2[1])) ||
      (max(P1[2], Q1[2]) < min(P2[2], Q2[2])) ||
      (min(P1[2], Q1[2]) > max(P2[2], Q2[2]))) {
    return(FALSE)
  }
  M <- cbind(Q1-P1, -Q2+P2)
  if (!(det(M) == 0)) {
    lambda <- solve(M) %*% (P2-P1)
    if (all((lambda<1) & (lambda>0))) {
      return(list(lambda=lambda, R=(1-lambda[1])*P1 + lambda[1]*Q1))
    }
  }
  return(FALSE)
}

## Remove intersections between adjacent segements in a circular path
remove.intersections <- function(P) {
  N <- nrow(P)
  for(i in 1:N) {
    fi <- find.intersection(P[i,],            P[mod1(i+1, N),],
                            P[mod1(i+2, N),], P[mod1(i+3, N),])
    if (is.list(fi)) {
      print(fi)
      P[mod1(i+1, N),] <- fi$R
      P <- P[-mod1(i+2, N),] 
      return(remove.intersections(P))
    }
  }
  return(P)
}
