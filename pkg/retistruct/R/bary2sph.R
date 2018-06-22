bary2sph <- function(Ib, T, P) {
  if (!is.matrix(Ib$p)) {
    Ib$p <- as.matrix(Ib$p, ncol=3)
  }
  if (!is.matrix(T)) {
    T <- as.matrix(T, ncol=3)
  }
  return(.Call("bary2sph", as.integer(Ib$idx), Ib$p, matrix(as.integer(T), ncol=3), P, PACKAGE="retistruct"))
}
