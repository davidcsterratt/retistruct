##' Constructor for AnnotatedOutline object.
##'
##' @title Constructor for AnnotatedOutline object
##' @param o \code{Outline} object
##' @return AnnotatedOutline object, with extra fields for tears
##' (\code{V0}, \code{VF} and \code{VB}), lattitude of rim \code{phi0}
##' and index of fixed point \code{i0}.
##' @author David Sterratt
AnnotatedOutline <- function(o){
  a <- o
  class(a) <- c("annotatedOutline", class(o))
  ## Trick to make V0, VB and VF "named numeric" of length 0
  a$V0 <- c(x=0)[0]
  a$VB <- c(x=0)[0]
  a$VF <- c(x=0)[0]
  a$phi0 <- 0
  a$i0 <- 1
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
##' \code{VF} and \code{VB}
##' @author David Sterratt
getTear <- function(o, tid) {
  return(with(o, c(V0=V0[tid], VB=VB[tid], VF=VF[tid])))
}

##' Compute the parent relationships for a potential set of tears on
##' an \code{AnnotatedOutline}. The function throws an error if tears
##' overlap.
##'
##' @title Compute the parent relationships for a set of tears
##' @param o \code{AnnotatedOutline} object
##' @param V0 Apices of tears
##' @param VB Backward vertices of tears
##' @param VF Forward vertices of tears
##' @return List
##' \item{\code{Rset}}{the set of points on the rim}
##' \item{\code{TFset}}{list containing indicies of points in each foward tear}
##' \item{\code{TBset}}{list containing indicies of points in each backward tear}
##' \item{\code{h}}{correspondence mapping}
##' \item{\code{hf}}{correspondence mapping in forward direction for
##'         points on boundary}
##' \item{\code{hb}}{correspondence mapping in backward direction for
##'         points on boundary}
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

  if (M > 0) {
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
##' @param a \code{AnnotatedOutline} object
##' @param pids Vector of three point IDs to be added
##' @return \code{AnnotatedOutline} object
##' @author David Sterratt
addTear <- function(a, pids) {
  M <- labelTearPoints(a, pids)
  V0 <- c(a$V0, M["V0"])
  VF <- c(a$VF, M["VF"])
  VB <- c(a$VB, M["VB"])
  ## This call will throw an error if tears are not valid
  suppressMessages(computeTearRelationships(a, V0, VB, VF))
  a$V0 <- V0
  a$VF <- VF
  a$VB <- VB
  a <- ensureFixedPointInRim(a)
  return(a)
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

setFixedPoint <- function(o, i0, name) {
  o$i0 <- i0
  names(o$i0) <- name
  o <- ensureFixedPointInRim(o)
  return(o)
}

##' Ensure that the fixed point \code{i0} is in the rim, not a tear.
##'
##' @title Ensure that the fixed point is in the rim, not a tear
##' @param o \code{annotatedOutline} object
##' @return o \code{annotatedOutline} object in which \code{i0} may
##' have been changed. 
##' @author David Sterratt
ensureFixedPointInRim <- function(o) {
  suppressMessages(t <- computeTearRelationships(o, o$V0, o$VB, o$VF))
  Rset <- t$Rset
  i0 <- o$i0
  if (!(i0 %in% Rset)) {
    o$i0 <- with(o, Rset[which.min(abs(Rset - i0))])
    if (!is.null(names(o$i0))) {
      warning(paste(names(o$i0)[1], "point has been moved to be in the rim"))
      names(o$i0) <- names(i0)
    } else {
      message("Fixed point has been moved to be in the rim")
    }
  }
  return(o)
}

plot.flat.annotatedOutline <- function(a, axt="n", ylim=NULL, ...) {
  NextMethod()
  args <- list(...)
  plot.markup <- is.null(args$markup) || args$markup

  if (plot.markup) {
    with(a, {
      if (length(V0) > 0) {
        points(P[VF,,drop=FALSE], col="red", pch="+")
        segments(P[V0,1], P[V0,2], P[VF,1], P[VF,2], col="red")
        points(P[VB,,drop=FALSE], col="orange", pch="+")
        segments(P[V0,1], P[V0,2], P[VB,1], P[VB,2], col="orange")
        points(P[V0,,drop=FALSE], col="cyan", pch="+")
        text(P[V0,,drop=FALSE]+100, labels=1:length(V0), col="cyan")
      }
      if (!is.null(names(i0))) {
        text(P[i0,1], P[i0,2], substr(names(i0)[1], 1, 1))
      }
    })
  }
}

## Stitch together tears in an outline
##
## Input arguments:
## P     - the coordinates of points in a mesh, including those the outline
## gf    - the forward pointer list
## gb    - the backward pointer list
## V0    - indicies of the apex of each tear
## VF    - indicies of the forward vertex of each tear
## VB    - indicies of the backward vertex of each tear
## i0    - the index of the landmark; this needs to be in the rim
##
## The function returns a list contatining:
## Rset  - the set of points on the rim
## i0    - the index of the landmark
## P     - a new set of meshpoints
## V0    - indicies of the apex of each tear
## VF    - indicies of the forward vertex of each tear
## VB    - indicies of the backward vertex of each tear
## TFset - list containing indicies of points in each foward tear
## TBset - list containing indicies of points in each backward tear
## gf    - new forward pointer list
## gb    - new backward pointer list
## h     - correspondence mapping
## hf    - correspondence mapping in forward direction for
##         points on boundary
## hb    - correspondence mapping in backward direction for
##         points on boundary
##
stitch.outline <- function(a) {
  
  r <- computeTearRelationships(a, a$V0, a$VB, a$VF)

  ## If not set, set the landmark marker index. Otherwise
  ## check it
  Rset <- r$Rset
  if (!(a$i0 %in% Rset)) {
    print(a$i0)
    print(Rset)
    stop("Fixed Point is not in rim")
  }

  P <- a$P
  V0 <- a$V0
  VF <- a$VF
  VB <- a$VB
  TFset <- r$TFset
  TBset <- r$TBset
  gf <- a$gf
  gb <- a$gb
  hf <- r$hf
  hb <- r$hb
  h <- r$h
  
  ## Insert points on the backward tears corresponding to points on
  ## the forward tears
  sF <-      stitch.insert.points(P, V0, VF, VB, TFset, TBset,
                                               gf, gb, hf, hb, h,
                                               "Forwards")

  ## Insert points on the forward tears corresponding to points on
  ## the backward tears
  sB <- with(sF,
             stitch.insert.points(P, V0, VB, VF, TBset, TFset,
                                  gb, gf, hb, hf, h,
                                  "Backwards"))
  ## Extract data from object
  P <- sB$P
  gf <- sB$gb
  gb <- sB$gf
  hf <- sB$hb
  hb <- sB$hf
  h <- sB$h

  ## Link up points on rim
  h[Rset] <- hf[Rset]
  
  ## Make sure that there are no chains of correspondences
  while (!all(h==h[h])) {
   h <- h[h]
  }

  s <- merge(list(Rset=Rset, i0=a$i0,
                  VF=VF, VB=VB, V0=V0,
                  TFset=TFset, TBset=TBset,
                  P=P, h=h, hf=hf, hb=hb,
                  gf=gf, gb=gb), a)
  class(s) <- c("stitchedOutline", class(a))
  return(s)
}

## Inner function responsible for inserting the points
stitch.insert.points <- function(P, V0, VF, VB, TFset, TBset, gf, gb, hf, hb, h,
                                 dir) {
  M <- length(V0)                       # Number of tears
  ## Iterate through tears to insert new points
  for (j in 1:M) {
    ## Compute the total path length along each side of the tear
    Sf <- path.length(V0[j], VF[j], gf, hf, P)
    Sb <- path.length(V0[j], VB[j], gb, hb, P)
    message(paste("Tear", j, ": Sf =", Sf, "; Sb =", Sb))

    ## For each point in the forward path, create one in the backwards
    ## path at the same fractional location
    message(paste("  ", dir, " path", sep=""))
    for (i in setdiff(TFset[[j]], c(V0[j], VF[j]))) {
      sf <- path.length(V0[j], i, gf, hf, P)
      ## If the point isn't at the apex, insert a point
      if (sf > 0) {
        message(paste("    i =", i,
                                "; sf/Sf =", sf/Sf,
                                "; sf =", sf))
        for (k in TBset[[j]]) {
          sb <- path.length(V0[j], k, gb, hb, P)
          message(paste("      k =", format(k, width=4),
                                  "; sb/Sb =", sb/Sb,
                                  "; sb =", sb))
          if (sb/Sb > sf/Sf) {
            break;
          }
          k0 <- k
          sb0 <- sb
        }

        ## If this point does not point to another, create a new point
        if ((hf[i] == i)) {
          f <- (sf/Sf*Sb-sb0)/(sb-sb0)
          message(paste("      Creating new point: f =", f))
          p <- (1-f) * P[k0,] + f * P[k,]

          ## Find the index of any row of P that matches p
          n <- anyDuplicated(rbind(P, p), fromLast=TRUE) 
          if (n == 0) {
            ## If the point p doesn't exist
            P <- rbind(P, p)
            ## Update forward and backward pointers
            n <- nrow(P)                    # Index of new point
            gb[n]     <- k
            gf[n]     <- gf[k]
            gb[gf[k]] <- n
            gf[k]     <- n

            ## Update correspondences
            hf[n] <- n
            hb[n] <- n
            h[i] <- n
            h[n] <- n
          } else {
            message(paste("      Point", n, "already exists"))
            h[i] <- n
            h[n] <- n
          }
        } else {
          ## If not creating a point, set the point to point to the forward pointer 
          h[i] <- hf[i]
        }
      }
    } 
  }
  return(list(P=P, hf=hf, hb=hb, gf=gf, gb=gb, h=h))
}
