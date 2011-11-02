##' Construct an outline object. This sanitises the input points
##' \code{P}, as described below.
##'
##' @title Outline constructor
##' @param P The points of the outline. The last point is not repeated.
##' @param scale The length of one unit of \code{P} in
##' micrometres. When images are present, this is the length of the
##' side of a pixel in the image.
##' @param im An image as a \code{raster} object
##' @return An \code{Outline} object containing the following:
##' \item{\code{P}}{A N-by-2 matrix of points of the \code{Outline} arranged in anticlockwise order}
##' \item{\code{gf}}{For each row of \code{P}, the index of \code{P} that is next in the outline travelling anticlockwise (forwards)}
##' \item{\code{gb}}{For each row of \code{P}, the index of \code{P} that is next in the outline travelling clockwise (backwards)}
##' \item{\code{im}}{The image as a \code{raster} object}
##' \item{\code{scale}}{The length of one unit of \code{P} in micrometres}
##' @author David Sterratt
Outline <- function(P, scale=NA, im=NULL) {
  o <- list()
  o$P <- P
  o$h <- 1:nrow(P)
  t <- triangulate.outline(o, n=NA)
  class(o) <- "outline"
  o$P <- t$P
  o$gf <- t$gf
  o$gb <- t$gb
  o$im <- im
  o$scale <- scale
  return(o)
}

##' Plot flat outline. If the optional argument \code{plot.image} is
##' \code{TRUE} the image (if it is present) is displayed behind the
##' outline. If the optional argument \code{scalebar} is present and
##' numeric, a scale bar of length \code{scalebar} mm is plotted. If
##' \code{scalebar} is \code{FALSE} the scale bar is supressed, and by
##' default a scale bar of 1mm is drawn.
##'
##' @title Flat plot of outline
##' @param o \code{Outline} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @method plot.flat outline
##' @author David Sterratt
plot.flat.outline <- function(o, axt="n", ylim=NULL, ...) {
  args <- list(...)
  plot.image <- is.null(args$image) || args$image
  scalebar <- ifelse(!is.null(args$scalebar),
                     args$scalebar, 1) # Default value of scalebar 1mm
  scalebar <- ifelse(is.numeric(scalebar) && !is.null(o$scale), scalebar, FALSE)

  with(o, {
    s <- which(!is.na(gb))                # source index
    d <- na.omit(gb)                      # destination index
    
    if (plot.image && !is.null(o$im)) {
      xs <- 1:ncol(im)
      ys <- 1:nrow(im)
      plot(NA, NA, xlim=c(0, max(xs)), ylim=c(0, max(ys)), asp=1,
           xaxt=axt, yaxt=axt, bty="n",
           xlab="", ylab="")
      ## rasterImage crashes on some systems, but not others.
      rasterImage(im, 0, 0, ncol(im), nrow(im))
    } else {
      xs <- P[s,1]
      ys <- P[s,2]
      suppressWarnings(plot(xs, ys, asp=1,
                            pch=".", xaxt=axt, yaxt=axt, xlab="", ylab="",
                            bty="n", ylim=ylim,  ...))
    }
    suppressWarnings(segments(P[s,1], P[s,2], P[d,1], P[d,2],
                              col=getOption("outline.col"), ...))

    ## Plot scalebar if required. scalebar is length in mm.
    if (scalebar && !is.na(scale)) {
      sby <- min(ys) - 0.02*(max(ys) - min(ys))
      sblen <- 1000*scalebar/(scale)
      lines(c(max(xs) - sblen, max(xs)),c(sby, sby), lwd=2)
    }
  })
}

##' Create a triangulation of the \code{Outline} object \code{o}.  The
##' minimum number of triangles in the triangulation is specified by
##' \code{n}.
##' 
##' @title Triangulate outline
##' @param o \code{\link{Outline}} object
##' @param n Minimum number of points in the triangulation
##' @param suppress.external.steiner If \code{TRUE} prevent the
##' addition of points in the outline. This happens to maintain
##' triangle quality.
##' @return A \code{triangulatedOutline} object containing the
##' following fields:
##' \item{\code{P}}{The set of new points, with the existing points at the start}
##' \item{\code{T}}{The triangulation}
##' \item{\code{Cu}}{Unique set of M connections, as M*2 matrix}
##' \item{\code{h}}{Correspondances mapping}
##' \item{\code{A}}{Array containing area of each triangle}
##' \item{\code{L}}{Length of each connection}
##' \item{\code{A.signed}}{Signed area of each triangle}
##' \item{\code{A.tot}}{Total area of outline}
##' \item{\code{gf}}{Forward pointers}
##' \item{\code{gb}}{Backward pointers}
##' \item{\code{S}}{Segments (from \code{\link{triangulate}})}
##' \item{\code{E}}{Edges (from \code{\link{triangulate}})}
##' \item{\code{EB}}{Edge boundaries (from \code{\link{triangulate}})}
##' @author David Sterratt
triangulate.outline <- function(o, n=200,
                                suppress.external.steiner=FALSE) {
  P <- o$P
  g <- o$gf
  h <- o$h

  if (is.null(h)) {
    h=1:nrow(P)
  }
  ## By default, segments are outline of points in order
  S <- cbind(1:nrow(P), c(2:nrow(P), 1))
  if (!is.null(g)) {
    S <- pointers2segments(g)
  }
  ## Make initial triangulation
  out <- triangulate(pslg(V=P, S=S), Y=TRUE, j=TRUE, Q=TRUE)

  ## It can be that there are crossovers in the segments. The
  ## triangulate() routine will reveal this as segments that are not
  ## on a boundary. We get rid of these segments by re-triangulating,
  ## only using boundary segments
  out <- triangulate(pslg(V=out$V, S=out$S[out$SB==1,]), Y=TRUE, j=TRUE, Q=TRUE)
  
  ## Sometimes a point exists which only belongs to one segment. The
  ## point to which it is connected, is itself connected by three
  ## segments. We want to get rid of these points, and the easiest way
  ## is to triangulate without the naughty points.
  i.bad <- which(table(out$S)==1)
  if (length(i.bad) > 0) {
    warning(paste("Bad points:", paste(i.bad, collapse=" ")))
    out <- triangulate(pslg(V=P[-i.bad,], S=S), Y=TRUE, j=TRUE, Q=TRUE)
  }

  ## Now determine the area
  A.tot <- sum(with(out, tri.area(V, T)))

  ## Produce refined triangulation
  P <- out$V
  S <- out$S
  if (!is.na(n)) {
    out <- triangulate(pslg(V=P, S=S), a=A.tot/n, q=20,
                       Y=suppress.external.steiner, j=TRUE,
                       Q=TRUE)
  }
  if (any(P != out$V[1:nrow(P),])) {
    stop("Points changed in triangulation")
  }
  P <- out$V
  T <- out$T

  ## Create pointers from segments

  ## To ensure the correct orientaion, we use the fact that the
  ## triangles are all anticlockwise in orinentation, and that the
  ## orientation of the first row of the segment matrix determines the
  ## orientation of all the other rows.

  ## We therefore find the triangle which contains the first segment
  S <- out$S
  T1 <- which(apply(T, 1, function(x) {all(S[1,] %in% x)}))

  ## Then find out which of the vertices in the triangle is not the
  ## one we need
  i <- which((T[T1,] %in% S[1,]) == FALSE)
  if (i == 3) S[1,] <- T[T1,c(1,2)]
  if (i == 2) S[1,] <- T[T1,c(3,1)]
  if (i == 1) S[1,] <- T[T1,c(2,3)]

  ## Now create the pointers from the segments
  gf <- segments2pointers(S)
  gb <- gf
  gb[na.omit(gf)] <- which(!is.na(gf))
  Rset <- na.omit(gf)
  
  ## Derive edge matrix from triangulation
  Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  ## If we are in the business of refining triangles (i.e. specifying
  ## n), remove lines which join non-ajancent parts of the outline
  if (!is.na(n)) {
    for (i in 1:nrow(Cu)) {
      C1 <- Cu[i,1]
      C2 <- Cu[i,2]
      if (all(Cu[i,] %in% Rset)) {
        if (!((C1 == gf[C2]) ||
              (C2 == gf[C1]))) {
          ## Find triangles containing the line
          ## segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="yellow")
          Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
          print(paste("Non-adjacent points in rim connected by line:", C1, C2))
          print(paste("In triangle:", Tind))
          ## Find points T1 & T2 in the two triangles which are not common
          ## with the edge
          T1 <- setdiff(T[Tind[1],], Cu[i,])
          T2 <- setdiff(T[Tind[2],], Cu[i,])
          print(paste("Other points in triangles:", T1, T2))
          ## Create a new point at the centroid of the four verticies
          ## C1, C2, T1, T2
          p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
          P <- rbind(P, p)
          n <- nrow(P)
          ## Remove the two old triangles, and create the four new ones
          T[Tind[1],] <- c(n, C1, T1)
          T[Tind[2],] <- c(n, C1, T2)
          T <- rbind(T,
                     c(n, C2, T1),
                     c(n, C2, T2))
        }
      }
    }

    ## Add the new points to the correspondances vector
    h <- c(h, (length(h)+1):nrow(P))

    ## Create the edge matrix from the triangulation
    Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
    Cu <- Unique(Cu, TRUE)
  }

  ## Swap orientation of triangles which have clockwise orientation
  A.signed <- tri.area.signed(P, T)
  T[A.signed<0,c(2,3)] <- T[A.signed<0,c(3,2)]
  A <- abs(A.signed)
  
  ## Find lengths of connections
  L <- vecnorm(P[Cu[,1],] - P[Cu[,2],])

  ## Check there are no zero-length lines
  if (any(L==0)) {
    print("WARNING: zero-length lines")
  }

  t <- merge(list(P=P, T=T, Cu=Cu, h=h,  A=A, L=L,
                  A.signed=A.signed, A.tot=A.tot,
                  gf=gf, gb=gb, S=out$S, E=out$E, EB=out$EB), o)
  class(t) <- addClass("triangulatedOutline", o)
  return(t)
}

##' Simplify an outline object by removing verticies bordering short
##' edges while not encroaching on any of the outline. At present,
##' this is done by finding concave vertices. It is safe to remove
##' these, at the expense of increasing the area a bit.
##'
##' @title Simplify an outline object by removing short edges
##' @param o \code{outline} object to simplify
##' @param min.frac.length the minumum length as a fraction of the
##' total length of the outline. 
##' @param plot whether to display plotting or not during simplification
##' @return Simlified \code{outline} object
##' @author David Sterratt
simplify.outline <- function(o, min.frac.length=0.001, plot=FALSE) {
  P <- o$P
  N <- nrow(P)                        # Number of vertices
  Q <- rbind(P, P[1,])                # Convenience variable
  v <- diff(Q)                         # Vector of each edge
  l <- vecnorm(v)                     # Length of each edge
  ## Compute outer products at each vertex
  e <- extprod3d(cbind(v[c(N, 1:(N-1)),], 0), cbind(v, 0))[,3]
  
  ## Find short edges
  S <- l/sum(l) < min.frac.length

  ## Find indicies of points that can be removed.
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
    message(paste("simplify.outline: Removing vertex", i.rem[1]))
    return(simplify.outline(list(P=P[-i.rem[1],], scale=o$scale, im=o$im)))
  } else {
    return(Outline(P, o$scale, o$im))
  }
}
