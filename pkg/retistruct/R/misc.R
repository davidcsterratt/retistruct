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

## Return standard names for lists containing any un-named elements
stdnames <- function(l) {
  nl <- names(l)
  if (length(l) >= 1) {
    ll <- paste("n", 1:length(l), sep="")
    ## Need to identify which elements are already named, i.e.
    ## elements which are not empty, and which contain one alphabetic
    ## character
    named <- which(nl!="" & grepl("[[:alpha:]]", nl))
    ll[named] <- nl[named]
  } else {
    ll <- nl
  }
  return(ll)
}

##' Return a new version of the list in which any un-named elements
##' have been given standardised names
##' 
##' @param l the list with un-named elements 
##' @return The list with standardised names
##' @author David Sterratt
name.list <- function(l) {
  if (is.list(l)) {
    names(l) <- stdnames(l)
  }
  return(l)
}

##' Parse dependencies
##' @param deps Text produced by, e.g., installed.packages()["packagename","Suggests"]
##' @return Table with package column, relationship column and version number
##' @author David Sterratt
parse.dependencies <- function(deps) {
  z <- unlist(strsplit(deps, ",\n", fixed = TRUE))
  pat <- "^([^\\([:space:]]+)[[:space:]]*\\(([^\\)]+)\\).*"
  deps <- cbind(sub(pat, "\\1", z), sub(pat, "\\2", z), NA)
  noversion <- deps[, 1] == deps[, 2]
  deps[noversion, 2] <- NA
  pat <- "[[:space:]]*([[<>=]+)[[:space:]]+(.*)"
  deps[!noversion, 2:3] <- c(sub(pat, "\\1", deps[!noversion, 2]),
                             sub(pat, "\\2", deps[!noversion, 2]))
  return(deps)
}

