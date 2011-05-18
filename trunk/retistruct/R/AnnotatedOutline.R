##' Constructor for AnnotatedOutline object.
##'
##' @title Constructor for AnnotatedOutline object
##' @param o \code{Outline} object
##' @return AnnotatedOutline object, with extra fields for tears
##' (\code{V0}, \code{VF} and \code{VB}) and lattitude of rim
##' (\code{phi0}).
##' @author David Sterratt
AnnotatedOutline <- function(o){
  a <- o
  class(a) <- c("annotatedOutline", class(o))
  a$V0 <- c()
  a$VB <- c()
  a$VF <- c()
  a$phi0 <- 0
  return(a)
}

##' Label a set of three unlabelled points supposed to refer to the
##' apex and vertcies of a cut and tear with the V0 (Apex), VF
##' (forward vertex) and VB (backward vertex) labels.
##'
##' @title Label three outline point indicies as apicies and vertices of tear
##' @param m the vector of three indicies
##' @param o Outline object
##' @return Vector of indicies labelled with V0, VF and VB
##' @author David Sterratt
labelTearPoints <- function(o, m) {
  with(o, {
    ## Each row of this matrix is a permutation of the markers
    p <- rbind(c(1, 2, 3),
               c(1, 3, 2),
               c(2, 1, 3),
               c(2, 3, 1),
               c(3, 1, 2),
               c(3, 2, 1))

    ## For each permuation of V0, VF, VB, measure the sum of length in
    ## the forwards direction from V0 to VF and in the backwards
    ## direction from V0 to VB. The permuation with the minimum distance
    ## is the correct one.
    tplmin <- Inf                      # The minimum path length
    h <- 1:nrow(P)                     # identity correspondence mapping
                                        # used for measuring distances
                                        # (this effectively ignores
                                        # sub-tears, but this doesn't
                                        # matter)
    for (i in 1:nrow(p)) {
      V0 <- m[p[i,1]]
      VF <- m[p[i,2]]
      VB <- m[p[i,3]]
      tpl <- path.length(V0, VF, gf, h, P) + path.length(V0, VB, gb, h, P)
      if (tpl < tplmin) {
        M <- m[p[i,]]
        tplmin <- tpl
      }
    }
    names(M) <- c("V0", "VF", "VB")
    return(M)})
}

##' Return index of tear in an AnnotatedOutline in which a point
##' appears
##'
##' @title Return index of tear
##' @param o \code{AnnotatedOutline} object
##' @param pid ID of point
##' @return ID of tear
##' @author David Sterratt
whichTear <- function(o, pid) {
  M <- with(o, cbind(V0, VF, VB))       # Tear matrix
  tid <- which(apply(pid==M, 1, any))[1]
  if (!length(tid))
    tid <- NA
  return(tid)
}

##' Return indicies of tear in AnnotatedOutline
##'
##' @title Return indicies of tear in AnnotatedOutline
##' @param o \code{AnnotatedOutline} object
##' @param tid Tear ID, which can be returned from \code{whichTear()}
##' @return Vector of three point IDs, labelled with \code{V0},
##' \code{VF} and \ccode{VB}
##' @author David Sterratt
getTear <- function(o, tid) {
  return(with(o, c(V0=V0[tid], VB=VB[tid], VF=VF[tid])))
}

##' Add tear to an AnnotatedOutline
##'
##' @title Add tear to an AnnotatedOutline
##' @param o \code{AnnotatedOutline} object
##' @param pids Vector of three point IDs to be added
##' @return \code{AnnotatedOutline} object
##' @author David Sterratt
addTear <- function(o, pids) {
  M <- labelTearPoints(r, pids)
  o$V0 <- c(o$V0, M["V0"])
  o$VF <- c(o$VF, M["VF"])
  o$VB <- c(o$VB, M["VB"])
  return(o)
}

##' Remove tear from an AnnotatedOutline
##'
##' @title Remove tear from an AnnotatedOutline
##' @param o \code{AnnotatedOutline} object
##' @param tid Tear ID, which can be returned from \code{whichTear()}
##' @return \code{AnnotatedOutline} object
##' @author David Sterratt
removeTear <- function(o, tid) {
  if (!is.na(tid)) {
    o$V0 <- o$V0[-tid]
    o$VF <- o$VF[-tid]
    o$VB <- o$VB[-tid]
  }
  return(o)
}

##' Given a tear matrix T with columns "V0", "VF", and "VB", check
##' that all tears are correct.
##'
##' @title Check that tears are all in the correct direction
##' @param o \code{AnnotatedOutline} object
##' @return If all is OK, returns empty vector.  If not, returns
##' indicies of problematic tears.
##' @author David Sterratt
checkTears <- function(o) {
  out <- c()
  if (length(o$V0)) {
    for (i in 1:length(o$V0)) {
      ## Extract the markers for this row
      m <- with(o, c(V0[i], VF[i], VB[i]))
      M <- labelTearPoints(o, m)
      if (!all(M == m)) {
        out <- c(out, i)
      }
    }
  }
  return(out)
}
