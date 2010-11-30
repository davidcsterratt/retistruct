misc.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
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

## Return standard names for lists containing any un-named elements
stdnames <- function(l) {
  nl <- names(l)
  if (length(l) >= 1) {
    ll <- paste("n", 1:length(l), sep="")
    named <- which(nl!="" & is.na(strtoi(nl)))
    ll[named] <- nl[named]
  } else {
    ll <- nl
  }
  return(ll)
}

## Return a new version of the list in which any un-named elements
## have been given standardised names
name.list <- function(l) {
  if (is.list(l)) {
    names(l) <- stdnames(l)
  }
  return(l)
}
