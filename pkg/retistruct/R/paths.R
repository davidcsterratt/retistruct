## Return next index in path
path.next <- function(i, g, h) {
  return(ifelse(h[i]==i, g[i], h[i]))
}

## Return sequence of indices in path between i and j, governed by
## pointer vector p
path <- function(i, j, g, h) {
  if (length(g) != length(h)) {
    stop("g and h must be the same length")
  }
  if (i == j) {
    return(i)
  } else {
    return(c(i, path(path.next(i, g, h), j, g, h)))
  }
}

## Return sequence of indices in path between i and j, governed by
## pointer vector p
path.length <- function(i, j, g, h, P) {
  if (any(is.na(c(i, j)))) {
    stop("i or j contains NA")
  }
  if (i == j) {
    return(0)
  } else {
    if (h[i] == i) {
      return(sqrt(sum((P[i,] - P[g[i],])^2)) + path.length(g[i], j, g, h, P))
    } else {
      return(path.length(h[i], j, g, h, P))
    }
  }
}
