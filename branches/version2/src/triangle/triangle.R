if (is.loaded("R_triangulate")) dyn.unload("R_triangle.so")
dyn.load("R_triangle.so")

triangulate <- function(P, a=NULL) {
  if (ncol(P) == 2) {
    P <- t(P)
  }
  PB <- rep(1, ncol(P))
  S <- rbind(1:ncol(P), c(2:ncol(P),1))
  SB <- rep(1, ncol(S))
  out <- .Call("R_triangulate",
               P,
               as.integer(PB),
               as.integer(S),
               as.integer(SB),
               a)
  names(out) <- c("Q", "B", "T")
  return(out)
}
