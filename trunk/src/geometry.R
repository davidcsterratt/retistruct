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
