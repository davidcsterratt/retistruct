if (is.loaded("R_triangulate")) dyn.unload("R_triangle.so")
dyn.load("R_triangle.so")

triangulate <- function(P, B) {
  out <- .Call("R_triangulate", P, as.integer(B))
  names(out) <- c("Q", "T")
  return(out)
}
