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

##' Return a new version of the list in which any unnamed elements
##' have been given standardised names
##' 
##' @param l the list with unnamed elements 
##' @return The list with standardised names
##' @author David Sterratt
name.list <- function(l) {
  if (is.list(l)) {
    names(l) <- stdnames(l)
  }
  return(l)
}

##' Parse dependencies
##' @param deps Text produced by, e.g., \code{installed.packages()["packagename","Suggests"]}
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

##' Interpolate values in image
##' @param im image to interpolate
##' @param P N by 2 matrix of x, y values at which to interpolate. x
##'   is in range \code{[0, ncol(im)]} and y is in range \code{[0, nrow(im)]}
##' @param invert.y If \code{FALSE} (the default), the y coordinate is
##'   zero at the top of the image. \code{TRUE} the zero y coordinate
##'   is at the bottom.
##' @return Vector of N interpolated values
##' @author David Sterratt
interpolate.image <- function(im, P, invert.y=FALSE) {
  N <- ncol(im)
  M <- nrow(im)
  x <- P[,1]
  y <- P[,2]
  if (invert.y) {
    y <- M - y
  }
  if (max(x) > N) {
    stop(paste0("max X-value of P (", max(x),") is bigger than number of cols in im (", N, ")"))
  }
  if (min(x) < 0) {
    stop(paste0("min X-value of P (", min(x),") is less than number of cols in im (", N, ")"))
  }
  if (max(y) > M) {
    stop(paste0("max Y-value of P (", max(y),") is bigger than number of rows in im (", M, ")"))
  }
  if (min(y) < 0) {
    stop(paste0("min Y-value of P (", min(y),") is less than number of rows in im (", M, ")"))
  }

  ## Assume that centres of pixels in image are at {(0.5, 0.5), (1.5,
  ## 0.5)...}  Find pixels whose centres surround each point in P. The
  ## four pixels will have indices in image of (i1, j1), (i1, j2),
  ## (i2, j1), (i2, j2). If a coordinate is less than 0.5, the index
  ## returned is 1. If it is more than the the maximum pixel index
  ## minus 0.5, then the index is M or N.
  j1 <- pmax(floor(x + 0.5), 1)
  j2 <- pmin(ceiling(x + 0.5), N)
  i1 <- pmax(floor(y + 0.5), 1)
  i2 <- pmin(ceiling(y + 0.5), M)

  ## Bilinear interpolation
  return(mapply(function(x, y, i1, i2, j1, j2) {
    cbind(j1 + 0.5 - x, x - j1 + 0.5) %*%
      t(im[c(i1,i2),c(j1,j2)]) %*%
      rbind(i1 + 0.5 - y, y - i1 + 0.5)
  }, x, y, i1, i2, j1, j2))
}

##' Reporting utility function
##'
##' Calls function speified by option \code{retistruct.report}
##' @param ... Arguments to reporting function
##' @author David Sterratt
##' @export
report <- function(...) {
  f <- getOption("retistruct.report")
  f(...)
}
