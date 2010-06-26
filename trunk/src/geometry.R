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
