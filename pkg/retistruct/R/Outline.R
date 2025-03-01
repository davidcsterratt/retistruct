##' Class containing basic information about flat outlines
##'
##' @description An Outline has contains the polygon describing the
##'   outline and an image associated with the outline.
##' @author David Sterratt
##' @importFrom R6 R6Class
Outline <- R6Class("Outline",
  inherit = OutlineCommon,
  public = list(
    ##' @field P A N-by-2 matrix of points of the \code{Outline}
    ##'   arranged in anticlockwise order
    P=NULL,
    ##' @field scale The length of one unit of \code{P} in the X-Y plane in arbitrary units
    scale=NULL,
    ##' @field scalez The length of one unit of \code{P} in the Z-direction in arbitrary units
    scalez=NULL,
    ##' @field units String giving units of scaled P, e.g. \dQuote{um}
    units=NA,
    ##' @field gf For each row of \code{P}, the index of \code{P} that
    ##'   is next in the outline travelling anticlockwise (forwards)
    gf=NULL,
    ##' @field gb For each row of \code{P}, the index of \code{P} that
    ##'   is next in the outline travelling clockwise (backwards)
    gb=NULL,
    ##' @field h For each row of \code{P}, the cut of that
    ##'   point (which will be to itself initially)
    h=NULL,
    ##' @field im An image as a \code{raster} object
    im=NULL,
    ##' @field dm Depthmap, with same dimensions as \code{im}, which indicates
    ##' height of each pixel in Z-direction
    dm=NULL,
    ##' @field A.fragments Areas of fragments
    A.fragments = NULL,

    ##' @field dm.inferna.window.min Minimum window size (in pixels) for inferring missing values in depthmaps
    dm.inferna.window.min = 10,
    ##' @field dm.inferna.window.max Minimum window size (in pixels) for inferring missing values in depthmaps
    dm.inferna.window.max = 100,
    ##' @description Construct an outline object. This sanitises the
    ##'   input points \code{P}.
    ##' @param fragments A list of N-by-2 matrix of points for each fragment of the \code{Outline}
    ##' @param scale The length of one unit of \code{P} in arbitrary units
    ##' @param im The image as a \code{raster} object
    ##' @param scalez The length of one unit of \code{P} in the Z-direction in arbitrary units. If \code{NA}, the depthmap is ignored.
    ##' @param dm Depthmap, with same dimensions as \code{im}, which indicates
    ##' height of each pixel in Z-direction
    ##' @param units String giving units of scaled P, e.g. \dQuote{um}
    initialize=function(fragments=list(), scale=NA, im=NULL, scalez=NA, dm=NULL, units=NA) {
      self$P <- matrix(0, 0, 4)
      colnames(self$P) <- c("X", "Y", "Z", "FID")
      if (!is.null(dm) & !is.null(im)) {
        if (all(dim(im) != dim(dm))) {
          stop("Image and depthmap must have the same dimensions")
        }
      }
      self$im <- im
      self$dm <- dm
      self$scale <- scale
      self$scalez <- scalez
      self$units <- sub(units, "um", "\U00B5m")
      if (!is.list(fragments)) {
        fragments <- list(fragments)
      }
      if (length(fragments) > 0) {
        for (i in 1:length(fragments)) {
          fragment <- Fragment$new()
          fragment$initializeFromPoints(fragments[[i]])
          pids <- self$addPoints(fragment$P, i)
          self$mapFragment(fragment, pids)
          self$A.fragments[i] <- fragment$A.tot
        }
      }
    },
    ##' @description Image accessor
    ##' @return An image as a \code{raster} object
    getImage = function() {
      return(self$im)
    },
    ##' @description Image setter
    ##' @param im An image as a \code{raster} object
    replaceImage = function(im) {
      if (all(dim(im) != dim(self$dm))) {
        stop("Image and depthmap must have the same dimensions")
      }
      self$im <- im
    },
    ##' @description Map the point IDs of a \link{Fragment} on the
    ##'   point IDs of this Outline
    ##' @param fragment \link{Fragment} to map
    ##' @param pids Point IDs in Outline of points in \link{Fragment}
    mapFragment = function(fragment, pids) {
      if (length(fragment$gf) != length(pids)) {
        stop("Number of fragment indices being mapped is not equal to number of pids supplied")
      }
      self$gf <- self$mapPids(fragment$gf, self$gf, pids)
      self$gb <- self$mapPids(fragment$gb, self$gb, pids)
    },
    ##' @description Map references to points
    ##' @param x References to point indices in source
    ##' @param y References to existing point indices in target
    ##' @param pids IDs of points in point register
    ##' @return New references to point indices in target
    mapPids = function(x, y, pids) {
      y[pids[which(is.na(x))]] <- NA
      nna <- which(!is.na(x))
      y[pids[nna]] <- pids[x[nna]]
      return(y)
    },
    ##' @description Get depth of points P
    ##' @param P matrix containing unscaled X-Y coordinates of points
    ##' @return Vector of unscaled Z coordinates of points P
    getDepth = function(P) {
      if (!is.null(self$dm)) {
        Z <- interpolate.image(self$dm, P, invert.y=TRUE,
                               wmin=self$dm.inferna.window.min,
                               wmax=self$dm.inferna.window.max)
      } else {
        Z <- rep(0, nrow(P))
      }
      return(Z)
    },

    ##' @description Add points to the outline register of points
    ##' @param P 2 or 3 column matrix of points to add
    ##' @param fid ID of fragment to which to add  the points
    ##' @return The ID of each added point in the register. If points already
    ##'   exist a point will not be created in the register,
    ##'   but an ID will be returned
    addPoints = function(P, fid) {
      if (!(is.vector(P) | is.matrix(P))) {
        stop("P must be a matrix with 2 or 3 columns, or a vector of length 2 or 3")
      }
      if (is.vector(P)) {
        if (length(P) == 2 | length(P) == 3) {
          P <- matrix(P, nrow=1)
        } else {
          stop("P should be vector of length 2 or 3")
        }
      }
      if (is.matrix(P)) {
        if (!(ncol(P) == 2 | ncol(P) == 3)) {
          stop("P should be a matrix with 2 or 3 columns")
        }
      }
      if (ncol(P) == 2) {
        P <- cbind(P, self$getDepth(P))
      }
      pids <- rep(NA, nrow(P))
      for (i in (1:nrow(P))) {
        if (nrow(self$P) == 0) {
          self$P <- rbind(self$P, c(P[i,], fid))
          pids[i] <- nrow(self$P)
          self$h[1] <- 1
        } else {
          ## Check point doesn't already exist
          id <- which(apply(t(self$P[,1:2,drop=FALSE]) == P[i,1:2], 2, all))
          if (length(id) > 1) {
            stop(paste("Point register has duplicates", self$P[id,], collapse=", "))
          }
          ## Point does exist
          if (length(id) == 1) {
            pids[i] <- id
          }
          ## Point doesn't exist
          if (length(id) == 0) {
            self$P <- rbind(self$P, c(P[i,], fid))
            pids[i] <- nrow(self$P)
            self$h <- c(self$h, pids[i])
          }
        }
      }
      return(pids)
    },
    ##' @description Get the point IDs in a fragment
    ##' @param fid fragment id of the points
    ##' @return Vector of point IDs, i.e. indices of the rows in
    ##'   the matrices returned by \code{getPoints} and
    ##'   \code{getPointsScaled}
    getFragmentPointIDs = function(fid) {
      return(which(self$P[,"FID"] == fid))
    },
    ##' @description Get the points in a fragment
    ##' @param fid fragment id of the points
    ##' @return Vector of points
    getFragmentPoints = function(fid) {
      return(self$P[self$getFragmentPointIDs(fid),c("X", "Y", "Z")])
    },
    ##' @description Get fragment
    ##' @param fid Fragment ID
    ##' @return The  \code{\link{Fragment}} object with ID fid
    getFragment = function(fid) {
      fpids <- self$getFragmentPointIDs(fid)
      map <- NULL
      map[fpids] <- 1:length(fpids)
      fragment <- Fragment$new()
      fragment$P = self$getFragmentPoints(fid)
      fragment$gf = map[self$gf[fpids]]
      fragment$gb = map[self$gb[fpids]]
      fragment$h = map[self$h[fpids]]
      return(fragment)
    },
    ##' @description Get fragment IDs from point IDS
    ##' @param pids Vector of point IDs
    ##' @return The Fragment ID to which each point belongs
    getFragmentIDsFromPointIDs = function(pids) {
      return(self$P[pids,"FID"])
    },
    ##' @description Get fragment IDs
    ##' @return IDs of all fragments
    getFragmentIDs = function() {
      return(unique(self$P[,"FID"]))
    },
    ##' @description Get unscaled mesh points
    ##' @param pids IDs of point to return
    ##' @return  Matrix with columns \code{X}, \code{Y} and \code{Z}
    getPoints = function(pids=NULL) {
      if (!is.null(pids)) {
        return(self$P[pids,c("X", "Y", "Z")])
      }
      return(self$P[,c("X", "Y", "Z")])
    },
    ##' @description Get X-Y coordinates of unscaled mesh points
    ##' @param pids IDs of point to return
    ##' @return  Matrix with columns \code{X} and \code{Y}
    getPointsXY = function(pids=NULL) {
      if (!is.null(pids)) {
        return(self$P[pids,c("X", "Y")])
      }
      return(self$P[,c("X", "Y")])
    },
    ##' @description Get scaled mesh points
    ##' @return  Matrix with columns \code{X} and \code{Y} which is
    ##'   exactly \code{scale} times the matrix returned by \code{getPoints}
    getPointsScaled = function() {
      if (is.na(self$scale)) {
        return(self$P[,c("X", "Y", "Z")])
      }
      if (is.na(self$scalez)) {
        return(cbind(self$scale*self$P[,c("X", "Y")],
                     Z=0))
      }
      return(cbind(self$scale*self$P[,c("X", "Y")],
                   Z=self$scalez*self$P[,"Z"]))
    },
    ##' @description Get set of points on rim
    ##' @return Vector of point IDs, i.e. indices of the rows in
    ##'   the matrices returned by \code{getPoints} and
    ##'   \code{getPointsScaled}
    getRimSet = function() {
      return(1:nrow(self$P))
    },
    ##' @description Get points on the edge of the outline
    ##' @return Vector of points IDs on outline
    getOutlineSet = function() {
      return(which(!is.na(self$gf)))
    },
    ##' @description Get lengths of edges of the outline
    ##' @return Vector of lengths of edges connecting neighbouring points
    getOutlineLengths = function() {
      return(vecnorm(self$getPointsScaled()[self$getOutlineSet(),] -
                     self$getPointsScaled()[self$gf[self$getOutlineSet()],]))
    },
    ##' @description Add a \link{FeatureSet}, e.g. a \link{PointSet}
    ##'   or \link{LandmarkSet}
    ##' @param fs \link{FeatureSet} to add
    addFeatureSet = function(fs) {
      if (fs$type %in% self$getFeatureSetTypes()) {
        stop(paste("There is already a", fs$type, "attached to this outline"))
      }
      self$featureSets <- c(self$featureSets, fs)
    }
  )
)

##' @method plot Outline
plot.Outline <- function(x, ...) {
  plot(x$self$P)
}

##' Plot flat \code{\link{Outline}}.
##'
##' @title Flat plot of outline
##' @param x \code{\link{Outline}} object
##' @param axt whether to plot axes
##' @param xlim x limits
##' @param ylim y limits
##' @param add If \code{TRUE}, don't draw axes; add to existing plot.
##' @param image If \code{TRUE} the image (if it is present) is
##' displayed behind the outline
##' @param scalebar If  numeric and if the Outline has a \code{scale}
##' field, a scale bar of length \code{scalebar} mm is plotted.  If
##' \code{scalebar} is \code{FALSE} or there is no scale information
##' in the \code{\link{Outline}} \code{x} the scale bar is suppressed.
##' @param rimset If \code{TRUE}, plot the points computed to be in the rim in the colour specified by the option \code{rimset.col}
##' @param pids If \code{TRUE}, plot point IDs
##' @param pid.joggle Amount to joggle point IDs by randomly
##' @param lwd.outline Line width of outline
##' @param ... Other plotting parameters
##' @method flatplot Outline
##' @author David Sterratt
##' @importFrom stats na.omit runif
##' @importFrom graphics lines plot segments points rasterImage
##' @export
flatplot.Outline <- function(x, axt="n",
                             xlim=NULL,
                             ylim=NULL,
                             add=FALSE,
                             image=TRUE,
                             scalebar=1,
                             rimset=FALSE,
                             pids=FALSE,
                             pid.joggle=0,
                             lwd.outline=1,
                             ...) {
  plot.image <- image
  ## If there is no scale information, don't present a scale bar
  scalebar <- ifelse(is.numeric(scalebar) && x$scale, scalebar, FALSE)


  s <- which(!is.na(x$gb))                # source index
  d <- na.omit(x$gb)                      # destination index

  if (!add) {
    im <- x$getImage()
    if (plot.image && !is.null(im)) {
      xs <- 1:ncol(im)
      ys <- 1:nrow(im)

      if (is.null(xlim)) {
        xlim <- c(0, max(xs))
      }
      if (is.null(ylim)) {
        ylim <- c(0, max(ys))
      }

      plot(NA, NA, xlim=xlim, ylim=ylim, asp=1,
           xaxt=axt, yaxt=axt, bty="n",
           xlab="", ylab="")
      ## rasterImage crashes on some systems, but not others.
      rasterImage(im, 0, 0, ncol(im), nrow(im))
    } else {
      xs <- x$getPoints()[s,"X"]
      ys <- x$getPoints()[s,"Y"]

      plot(xs, ys, asp=1,
           pch=".", xaxt=axt, yaxt=axt, xlab="", ylab="",
           bty="n", xlim=xlim, ylim=ylim)
    }
  }
  if (rimset) {
    points(x$P[x$getRimSet(),], col=getOption("rimset.col"), pch=19)
  }
  segments(x$P[s,1], x$P[s,2], x$P[d,1], x$P[d,2],
           col=getOption("outline.col"), lwd=lwd.outline)
  if (pids) {
    text(x$P[,"X"], x$P[,"Y"] + runif(nrow(x$P), -pid.joggle, pid.joggle), 1:nrow(x$P), ...)
  }

  ## Plot scalebar if required. scalebar is length in mm.
  if (!add && scalebar && !is.na(x$scale)) {
    sby <- min(ys) - 0.02*(max(ys) - min(ys))
    sblen <- 1000*scalebar/(x$scale)
    lines(c(max(xs) - sblen, max(xs)),c(sby, sby), lwd=2)
  }
}

##' Simplify a outline object by removing vertices bordering short
##' edges while not encroaching on any of the outline. At present,
##' this is done by finding concave vertices. It is safe to remove
##' these, at the expense of increasing the area a bit.
##'
##' @title Simplify an outline object by removing short edges
##' @param P points to simplify
##' @param min.frac.length the minimum length as a fraction of the
##' total length of the outline.
##' @param plot whether to display plotting or not during simplification
##' @return Simplified \code{outline} object
##' @author David Sterratt
simplifyOutline <- function(P, min.frac.length=0.001, plot=FALSE) {
  N <- nrow(P)                        # Number of vertices
  Q <- rbind(P, P[1,])                # Convenience variable
  v <- diff(Q)                         # Vector of each edge
  l <- vecnorm(v)                     # Length of each edge
  ## Compute outer products at each vertex
  e <- extprod3d(cbind(v[c(N, 1:(N-1)),], 0), cbind(v, 0))[,3]

  ## Find short edges
  S <- l/sum(l) < min.frac.length

  ## Find indices of points that can be removed.
  ## They have to be concave or colinear (e<=0). And they need to border short edges
  i.rem <- which((e <= 0) & (S | (S[c(N, 1:(N-1))])))

  if (plot) {
    ## Plot short edges...
    plot(P, col="white")
    if (any(S)) {
      segments(P[S,1], P[S,2], P[S,1]+v[S,1], P[S,2] + v[S,2], col="red")
    }
    ## and longer ones.
    if (any(!S)) {
      segments(P[!S,1], P[!S,2], P[!S,1]+v[!S,1], P[!S,2] + v[!S,2], col="black")
    }
    ## Plot colinear, convex and concave points
    points(P[e>0,1], P[e>0, 2], col="green")
    points(P[e==0,1], P[e==0, 2], col="orange")
    points(P[e<0,1], P[e<0, 2], col="blue")
    ## Plot points to remove
    points(P[i.rem,1], P[i.rem, 2], pch="X", col="brown")
  }

  ## If there are any points to remove, remove the first one
  if (length(i.rem) > 0) {
    message(paste("simplifyOutline: Removing vertex", i.rem[1]))
    return(simplifyOutline(P[-i.rem[1],],
                            min.frac.length=min.frac.length, plot=plot))
  } else {
    return(P)
  }
}
