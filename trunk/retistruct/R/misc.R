misc.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev: 277$"))))
}

## mod1(i, N)
##
## A modulus function that returns numbers in the
## range 1 to N, rather than 0 to N-1
##
## Specifically, it returns
## i , if i <= N
## i - N, if i > N
## 
mod1 <- function(i, N) {
  return((i - 1) %% N + 1)
}

## flipud(M)
##
## Flip the rows of matrix M, so that the bottom row becomes the
## top row
##
flipud <- function(M) {
  return(M[nrow(M):1,])
}

## merge.lists(l, m)
##
## Merge the contents of list l and list m, overwriting any names in
## l which also occur in m
##
merge.lists <- function(l, m) {
  for (n in names(m)) {
    l[n] = m[n]
  }
  return(l)
}
