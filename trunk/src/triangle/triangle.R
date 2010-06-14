dyn.load("R_triangle.so")
conv <- function(a, b)
       .C("convolve",
          as.double(a),
          as.integer(length(a)),
          as.double(b),
          as.integer(length(b)),
          ab = double(length(a) + length(b) - 1))$ab

conv2 <- function(a, b)
  .Call("convolve2", a, b)

triangulate <- function(P) {
  out <- .Call("R_triangulate", P)
  names(out) <- c("Q", "T")
  return(out)
}
