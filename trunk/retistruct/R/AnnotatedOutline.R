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

##' Compute the parent relationships for a potential set of tears on
##' an \code{AnnotatedOutline}. The function throws an error if tears
##' overlap
##'
##' @title Compute the parent relationships for a set of tears
##' @param o \code{AnnotatedOutline} object
##' @param V0 Apices of tears
##' @param VB Backward vertices of tears
##' @param VF Forward vertices of tears
##' @return List
##' \item{\code{Rset}} the set of points on the rim
##' \item{\code{TFset}} list containing indicies of points in each foward tear
##' \item{\code{TBset}} list containing indicies of points in each backward tear
##' \item{\code{h}} correspondence mapping
##' \item{\code{hf}} correspondence mapping in forward direction for
##'         points on boundary
##' \item{\code{hb}} correspondence mapping in backward direction for
##'         points on boundary
##' @author David Sterratt
computeTearRelationships <- function(o, V0, VB, VF) {
  ## Initialise the set of points in the rim
  ## We don't assume that P is the entire set of points; instead
  ## get this information from the pointer list.
  N <- nrow(o$P)                          # Number of points
  h <- 1:N                              # Initial correspondences
  hf <- h
  hb <- h
  M <- length(V0)                       # Number of tears
  i.parent <- rep(0, M)                 # Index of parent tear.
                                        # Is 0 if root otherwise
                                        # index of tear if in forward side
                                        # or negative index if in backward side 
  Rset <- na.omit(o$gf)
  
  ## Create lists of forward and backward tears
  TFset <- list()
  TBset <- list()
  
  ## Iterate through the tears to create tear sets and rim set
  for (j in 1:M) {
    ## Create sets of points for each tear and remove these points from
    ## the rim set
    ## message(paste("Forward tear", j))
    TFset[[j]] <- mod1(path(V0[j], VF[j], o$gf, h), N)
    TBset[[j]] <- mod1(path(V0[j], VB[j], o$gb, h), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))
  }
  
  ## Search for parent tears
  ## Go through all tears
  for (j in 1:M) {
    for (k in setdiff(1:M, j)) {
      ## If this tear is contained in a forward tear
      if (all(c(V0[j], VF[j], VB[j]) %in% TFset[[k]])) {
        i.parent[j] <- k
        message(paste("Tear", j, "child of forward side of tear", k))
        ## Set the forward pointer
        hf[VB[j]] <- VF[j]
        ## Remove the child tear points from the parent
        TFset[[k]] <- setdiff(TFset[[k]],
                              setdiff(c(TBset[[j]], TFset[[j]]), c(VB[j], VF[j])))
        ## message(TFset[[k]])
      } else {
        ## If this tear is contained in a backward tear
        if (all(c(V0[j], VF[j], VB[j]) %in% TBset[[k]])) {
          i.parent[j] <- -k
          message(paste("Tear", j, "child of backward side of tear", k))
          ## Set the forward pointer
          hb[VF[j]] <- VB[j]
          ## Remove the child tear points from the parent
          TBset[[k]] <- setdiff(TBset[[k]],
                                setdiff(c(TBset[[j]], TFset[[j]]), c(VB[j], VF[j])))
        } else {
          if (any(c(V0[j], VF[j], VB[j]) %in%
                  setdiff(union(TFset[[k]], TBset[[k]]), c(VF[k], VB[k])))) {
            stop(paste("Tear", j, "overlaps with tear", k))
          }
        }
      }
    }
    if (i.parent[j] == 0) {
      message(paste("Tear", j, "child of rim"))
      hf[VB[j]] <- VF[j]
      hb[VF[j]] <- VB[j]
    }
  }
  return(list(Rset=Rset,
              TFset=TFset,
              TBset=TBset,
              h=h,
              hf=hf,
              hb=hb))
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
  V0 <- c(o$V0, M["V0"])
  VF <- c(o$VF, M["VF"])
  VB <- c(o$VB, M["VB"])
  ## This call will throw an error if tears are not valid
  computeTearRelationships(o, V0, VB, VF)
  o$V0 <- V0
  o$VF <- VF
  o$VB <- VB
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

