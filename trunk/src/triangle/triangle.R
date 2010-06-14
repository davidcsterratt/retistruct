dyn.load("R_triangle.so")

triangulate <- function(P) {
  out <- .Call("R_triangulate", P)
  names(out) <- c("Q", "T")
  return(out)
}
