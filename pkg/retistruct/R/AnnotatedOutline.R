##' Class containing functions and data relating to annotating outlines
##'
##' @description An AnnotatedOutline contains a function to annotate 
##'   tears on the outline.
##'
##' @return AnnotatedOutline object, with extra fields for tears
##'   latitude of rim \code{phi0} and index of fixed point \code{i0}.
##' @author David Sterratt
##' @export
##' @importFrom R6 R6Class
##' @examples
##' P <- rbind(c(1,1),   c(2,1),  c(2,-1),
##'            c(1,-1),  c(1,-2), c(-1,-2),
##'            c(-1,-1), c(-2,-1),c(-2,1),
##'            c(-1,1),  c(-1,2), c(1,2))
##' o <- TriangulatedOutline$new(P)
##' o$addTear(c(3, 4, 5))
##' o$addTear(c(6, 7, 8))
##' o$addTear(c(9, 10, 11))
##' o$addTear(c(12, 1, 2))
##' flatplot(o)
AnnotatedOutline <- R6Class("AnnotatedOutline",
  inherit = PathOutline,
  public = list(
    ##' @field tears Matrix in which each row represents a tear by the
    ##'   indices into the outline points of the apex (\code{V0}) and
    ##'   backward (\code{VB}) and forward (\code{VF}) points
    tears = NULL,
    ##' @field phi0 rim angle in radians
    phi0 = 0,
    ##' @field lambda0 longitude of fixed point
    lambda0 = 0,
    ##' @field i0 index of fixed point
    i0 = 1,
    ##' @description Constructor
    ##' @param ... Parameters to \code{\link{PathOutline}}
    initialize = function(...) {
      super$initialize(...)
      self$tears <- matrix(0, 0, 3)
      colnames(self$tears) <- c("V0","VB","VF")
    },
    ##' @description Label a set of three unlabelled points supposed
    ##'   to refer to the apex and vertices of a cut and tear with the \code{V0}
    ##'   (Apex), \code{VF} (forward vertex) and \code{VB} (backward vertex) labels.
    ##' @param pids the vector of three indices
    ##' @return Vector of indices labelled with \code{V0}, \code{VF} and \code{VB}
    labelTearPoints = function(pids) {
      if (length(unique(pids)) != 3) {
        stop("Tears have to be defined by 3 points")
      }

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
      tplmin <- Inf                # The minimum path length
      h <- 1:nrow(self$P)          # identity correspondence mapping
                                        # used for measuring distances
                                        # (this effectively ignores
                                        # sub-tears, but this doesn't
                                        # matter)
      for (i in 1:nrow(p)) {
        V0 <- pids[p[i,1]]
        VB <- pids[p[i,2]]
        VF <- pids[p[i,3]]
        tpl <- path.length(V0, VF, self$gf, h, self$P) +
               path.length(V0, VB, self$gb, h, self$P)
        if (tpl < tplmin) {
          tear <- pids[p[i,]]
          tplmin <- tpl
        }
      }
      names(tear) <- c("V0", "VB", "VF")
      return(tear)

    },
    ##' @description Return index of tear in an AnnotatedOutline in which a point appears
    ##' @param pid ID of point
    ##' @return ID of tear
    whichTear = function(pid) {
      tid <- which(apply(pid==self$tears, 1, any))[1]
      if (!length(tid))
        tid <- NA
      return(tid)
    },
    ##' @description Return indices of tear in AnnotatedOutline
    ##' @param tid Tear ID, which can be returned from \code{whichTear()}
    ##' @return Vector of three point IDs, labelled with \code{V0},
    ##' \code{VF} and \code{VB}
    getTear = function(tid) {
      if (tid > nrow(self$tears)) {
        return(NA)
      }
      return(c(self$tears[tid,]))
    },
    ##' @description Get tears
    ##' @return Matrix of tears
    getTears = function() {
      return(self$tears)
    },
    ##' @description Compute the parent relationships for a potential
    ##'   set of tears. The function throws an error if tears overlap.
    ##' @param tears Matrix containing columns \code{V0} (Apices of tears) 
    ##'   \code{VB} (Backward vertices of tears) and \code{VF} (Forward
    ##'   vertices of tears)
    ##' @return List containing
    ##' \itemize{
    ##' \item{\code{Rset}}{the set of points on the rim}
    ##' \item{\code{TFset}}{list containing indices of points in each forward tear}
    ##' \item{\code{TBset}}{list containing indices of points in each backward tear}
    ##' \item{\code{h}}{correspondence mapping}
    ##' \item{\code{hf}}{correspondence mapping in forward direction for
    ##'         points on boundary}
    ##' \item{\code{hb}}{correspondence mapping in backward direction for
    ##'         points on boundary}
    ##' }
    computeTearRelationships = function(tears=NULL) {
      if (is.null(tears)) {
        tears <- self$getTears()
      }
      ## Initialise the set of points in the rim
      ## We don't assume that P is the entire set of points; instead
      ## get this information from the pointer list.
      N <- nrow(self$P)                 # number of points

      h <- 1:N                          # Initial correspondences
      hf <- h
      hb <- h
      M <- nrow(tears)        # Number of tears
      i.parent <- rep(0, M)   # Index of parent tear.
                                        # Is 0 if root otherwise
                                        # index of tear if in forward side
                                        # or negative index if in backward side 
      Rset <- na.omit(self$gf)
      
      ## Create lists of forward and backward tears
      TFset <- list()
      TBset <- list()
      ## Convenience variables
      V0 <- tears[,"V0"]
      VF <- tears[,"VF"]
      VB <- tears[,"VB"] 
      
      if (M > 0) {
        ## Iterate through the tears to create tear sets and rim set
        for (j in 1:M) {
          ## Create sets of points for each tear and remove these points from
          ## the rim set
          ## message(paste("Forward tear", j))
          TFset[[j]] <- mod1(path(V0[j], VF[j], self$gf, h), N)
          TBset[[j]] <- mod1(path(V0[j], VB[j], self$gb, h), N)
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
              report(paste("Tear", j, "child of forward side of tear", k))
              ## Set the forward pointer
              hf[VB[j]] <- VF[j]
              ## Remove the child tear points from the parent
              TFset[[k]] <- setdiff(TFset[[k]],
                                    setdiff(c(TBset[[j]], TFset[[j]]), c(VB[j], VF[j])))
              ## report(TFset[[k]])
            } else {
              ## If this tear is contained in a backward tear
              if (all(c(V0[j], VF[j], VB[j]) %in% TBset[[k]])) {
                i.parent[j] <- -k
                report(paste("Tear", j, "child of backward side of tear", k))
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
            report(paste("Tear", j, "child of rim"))
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
    },
    ##' @description Add tear to an AnnotatedOutline 
    ##' @param pids Vector of three point IDs to be added
    addTear = function(pids) {
      M <- self$labelTearPoints(pids)
      tears <- rbind(self$tears,  M[c("V0","VB","VF")])
      ## This call will throw an error if tears are not valid
      suppressMessages(self$computeTearRelationships(tears))
      self$tears <- tears
      self$ensureFixedPointInRim()
    },
    ##' @description Remove tear from an AnnotatedOutline
    ##' @param tid Tear ID, which can be returned from \code{whichTear()}
    removeTear = function(tid) {
      if (!is.na(tid)) {
        self$tears <- self$tears[-tid,]
      }
    },
    ##' @description Check that all tears are correct.
    ##' @return If all is OK, returns empty vector.  If not, returns
    ##'   indices of problematic tears.
    checkTears = function() {
      out <- c()
      if (nrow(self$tears)) {
        for (i in 1:nrow(self$tears)) {
          ## Extract the markers for this row
          if (!all(self$getTear(i) == self$labelTearPoints(self$getTear(i)))) {
            out <- c(out, i)
          }
        }
      }
      return(out)
    },
    ##' @description Set fixed point
    ##' @param i0 Index of fixed point
    ##' @param name Name of fixed point
    setFixedPoint = function(i0, name) {
      self$i0 <- i0
      names(self$i0) <- name
      self$ensureFixedPointInRim()
    },
    ##' @description Get point ID of fixed point
    ##' @return Point ID of fixed point
    getFixedPoint = function() {
      return(self$i0)
    },
    ##' @description Get point IDs of points on rim
    ##' @return Point IDs of points on rim
    getRimSet = function() {
      TR <- self$computeTearRelationships(self$tears)
      return(TR$Rset)
    },
    ##' @description Ensure that the fixed point \code{i0} is in the rim, not a tear.
    ##' Alters object in which \code{i0} may have been changed. 
    ensureFixedPointInRim = function() {
      Rset <- self$getRimSet()
      i0 <- self$i0
      if (!(i0 %in% Rset)) {
        self$i0 <- Rset[which.min(abs(Rset - self$i0))]
        if (!is.null(names(self$i0))) {
          warning(paste(names(self$i0)[1], "point has been moved to be in the rim"))
          names(self$i0) <- names(i0)
        } else {
          report("Fixed point has been moved to be in the rim")
        }
      }
    },
    ##' @description Get lengths of edges on rim
    ##' @return Vector of rim lengths
    getRimLengths = function() {
      ## Rim set
      rs <- self$getRimSet()
      ## Destination points
      d <- self$nextPoint(rs)
      ## Source and destination points have to be in rim set
      s <- rs[d %in% rs]
      d <- d[d %in% rs]
      return(vecnorm(self$getPointsScaled()[d,] -
                     self$getPointsScaled()[s,]))
    }
  )
)

##' Plot flat \code{\link{AnnotatedOutline}}. The user markup is
##' displayed by default. 
##'
##' @title Flat plot of AnnotatedOutline
##' @param x \code{\link{AnnotatedOutline}} object
##' @param axt whether to plot axes
##' @param xlim x-limits
##' @param ylim y-limits
##' @param markup If \code{TRUE}, plot markup
##' @param ... Other plotting parameters
##' @method flatplot AnnotatedOutline
##' @importFrom graphics text
##' @author David Sterratt
##' @export
flatplot.AnnotatedOutline <- function(x, axt="n",
                                      xlim=NULL,
                                      ylim=NULL,
                                      markup=TRUE,
                                      ...) {
  NextMethod()

  if (markup) {
    tears <- x$getTears()
    P <- x$P
    i0 <- x$i0
    gf <- x$gf
    h <- 1:nrow(x$P)
    if (nrow(tears) > 0) {
      V0 <- tears[,"V0"]
      VB <- tears[,"VB"]
      VF <- tears[,"VF"]
      points(P[VF,,drop=FALSE], col=getOption("TF.col"), pch="+")
      segments(P[V0,1], P[V0,2], P[VF,1], P[VF,2], col=getOption("TF.col"))
      points(P[VB,,drop=FALSE], col=getOption("TB.col"), pch="+")
      segments(P[V0,1], P[V0,2], P[VB,1], P[VB,2], col=getOption("TB.col"))
      points(P[V0,,drop=FALSE], col=getOption("V.col"), pch="+")
      text(P[V0,,drop=FALSE] + 0.02*(max(P[,1])-min(P[,1])),
           labels=1:length(V0), col=getOption("V.col"))
      segments(P[VF,1], P[VF,2], P[VB,1], P[VB,2], col=getOption("TF.col"), lty=2)
    }
    if (!is.null(names(i0))) {
      text(P[i0,1], P[i0,2], substr(names(i0)[1], 1, 1), col=getOption("V.col"))
    }
  }
}
