##' Reconstruct outline into spherical surface. Reconstruction
##' proceeds in a number of stages:
##'
##' \enumerate{
##' 
##' \item The flat object is triangulated with at least \code{n}
##' triangles. This can introduce new vertices in the rim. 
##'
##' \item The triangulated object is stitched.
##'
##' \item The stitched object is triangulated again, but this time it
##' is not permitted to add extra vertices to the rim.
##'
##' \item The corresponding points determined by the stitching process
##' are merged to form a new set of merged points and a new
##' triangulation.
##'
##' \item The merged points are projected roughly to a sphere.
##'
##' \item The locations of the points on the sphere are moved so as to
##' minimise the energy function.
##' }
##'
##' @title Reconstruct outline into spherical surface
##' @importFrom geometry tsearch sph2cart
##' @author David Sterratt
ReconstructedOutline <- R6Class("ReconstructedOutline",
  inherit = OutlineCommon,
  public = list(
    ol = NULL,                            # Annotated outline
    Pt = NULL,
    Tt = NULL,
    Ct = NULL,
    Cut = NULL,
    Bt = NULL,
    Lt = NULL,
    ht = NULL,
    u = NULL,
    U = NULL,
    Rsett = NULL,
    i0t = NULL,
    H = NULL,
    Ht = NULL,
    phi0 = NULL,
    R = NULL,
    lambda0 = NULL,
    lambda = NULL,
    phi = NULL,
    Ps = NULL,
    n = 500,
    alpha = 8,
    x0 = 0.5,
    nflip0 = NULL,
    nflip = NULL,
    opt = NULL,
    E.tot = NULL,
    E.l = NULL,
    mean.strain = NULL,
    mean.logstrain = NULL,
    ims = NULL,
    immask = NULL,
    report = NULL,
    debug = NULL,
    ## Optional function to transform FeatureSets linked to this outline
    featureSetTransform = function(fs) { return(fs) },
    ## @param o \code{\link{AnnotatedOutline}} object, containing the following information:\describe{
    ## \item{\code{P}}{outline points as N-by-2 matrix}
    ## \item{\code{V0}}{indices of the apex of each tear}
    ## \item{\code{VF}}{indices of the forward vertex of each tear}
    ## \item{\code{VB}}{indices of the backward vertex of each tear}
    ## \item{\code{i0}}{index of the landmark on the rim}
    ## \item{\code{phi0}}{latitude of rim of partial sphere}
    ## \item{\code{lambda0}}{longitude of landmark on rim}
    ## }
    ## @param n Number of points in triangulation.
    ## @param alpha Area scaling coefficient
    ## @param x0 Area cut-off coefficient
    ## @param plot.3d Whether to show 3D picture during optimisation.
    ## @param dev.flat Device to plot grid onto. Value of \code{NA} (default)
    ## means no plotting.
    ## @param dev.polar Device display projection. Value of NA
    ## (default) means no plotting.
    ## @param report Function to report progress.
    ## @param debug If \code{TRUE} print extra debugging output
    ## @return \code{reconstructedOutline} object containing the input
    ## information and the following modified and extra information:
    ## \item{\code{P}}{New set of points in flattened object}
    ## \item{\code{gf}}{New set of forward pointers in flattened object}
    ## \item{\code{gb}}{New set of backward pointers in flattened object}
    ## \item{\code{phi}}{latitude of new points on sphere}
    ## \item{\code{lambda}}{longitude of new points on sphere}
    ## \item{\code{Tt}}{New triangulation}
    initialize = function(ol,
                          n=500, alpha=8, x0=0.5,
                          plot.3d=FALSE, dev.flat=NA, dev.polar=NA, report=NULL,
                          debug=FALSE) {
      self$n <- n
      self$alpha <- alpha
      self$x0 <- x0
      if (is.null(report)) {
        self$report <- ol$report
      } else {
        self$report <- report
      }
      self$debug <- debug
      ol$triangulate()
      ol$stitchTears()
      ol$triangulate(suppress.external.steiner=TRUE)
      if (length(ol$corrs)) {
        ol$triangulate(suppress.external.steiner=TRUE)
      }
      ## Transform the rim set
      ## ol$orderRset()
      self$ol <- ol
      self$phi0 <- ol$phi0
      self$lambda0 <- ol$lambda0

      self$report("Merging points...")
      self$mergePointsEdges()

      self$report("Projecting to sphere...")
      self$projectToSphere()
    },
    reconstruct = function(plot.3d=FALSE, dev.flat=NA, dev.polar=NA) {

      ##   ## Initial plot in 3D space
      ##   if (plot.3d) {
      ##     sphericalplot(r)
      ##   }
      ## }

      ## Check for flipped triangles and record initial number
      ft <- flipped.triangles(self$phi, self$lambda, self$Tt, self$R)
      self$nflip0 <- sum(ft$flipped)

      self$report("Optimising mapping with no area constraint using BFGS...")
      self$optimiseMapping(alpha=0, x0=0, nu=1,
                           plot.3d=plot.3d,
                           dev.flat=dev.flat, dev.polar=dev.polar)
      self$report("Optimising mapping with area constraint using FIRE...")
      ## FIXME: Need to put in some better heuristics for scaling
      ## maxmove, and perhaps other parameters
      self$optimiseMappingCart(alpha=self$alpha, x0=self$x0, nu=1,
                               dtmax=500, maxmove=0.002*sqrt(self$ol$A.tot),
                               tol=1e-5,
                               plot.3d=plot.3d,
                               dev.flat=dev.flat, dev.polar=dev.polar)
      self$report("Optimising mapping with strong area constraint using BFGS...")
      self$optimiseMapping(alpha=self$alpha, x0=self$x0, nu=1,
                           plot.3d=plot.3d,
                           dev.flat=dev.flat, dev.polar=dev.polar)
      self$report("Optimising mapping with weak area constraint using BFGS...")
      self$optimiseMapping(alpha=self$alpha, x0=self$x0, nu=0.5,
                           plot.3d=plot.3d,
                           dev.flat=dev.flat, dev.polar=dev.polar)
      
      self$report(paste("Mapping optimised. Deformation energy E:", format(self$opt$value, 5),
                        ";", self$nflip, "flipped triangles."))
    },
    ## This function creates merged and transformed versions (all
    ## suffixed with \code{t}) of a number of existing variables, as well
    ## as a matrix \code{Bt}, which maps a binary vector representation
    ## of edge indices onto a binary vector representation of the
    ## indices of the points linked by the edge.
    ## @title  Merge stitched points and edges 
    ## @param t A \code{StitchedOutline} object in which points that have
    ## been added by stitching have been triangulated
    ## @return Adds following fields to input
    ## \item{\code{Pt}}{Transformed point locations}
    ## \item{\code{Tt}}{Transformed triangulation}
    ## \item{\code{Ct}}{Transformed connection set}
    ## \item{\code{Cut}}{Transformed symmetric connection set}
    ## \item{\code{Bt}}{Transformed binary vector representation
    ## of edge indices onto a binary vector representation of the
    ## indices of the points linked by the edge}
    ## \item{\code{Lt}}{Transformed edge lengths}
    ## \item{\code{ht}}{Transformed correspondences}
    ## \item{\code{u}}{Indices of unique points in untransformed space}
    ## \item{\code{U}}{Transformed indices of unique points in untransformed space}
    ## \item{\code{Rset}}{The set of points on the rim (which has been reordered)}
    ## \item{\code{Rsett}}{Transformed set of points on rim}
    ## \item{\code{i0t}}{Transformed index of the landmark}
    ## \item{H}{mapping from edges onto corresponding edges}
    ## \item{Ht}{Transformed mapping from edges onto corresponding edges}
    ## @author David Sterratt
    ## @export
    mergePointsEdges = function() {
      h <- self$ol$h
      T <- self$ol$T
      Cu <- self$ol$Cu
      L <- self$ol$L
      P <- self$ol$getPointsScaled()
      gf <- self$ol$gf
      
      ## Form the mapping from a new set of consecutive indices
      ## the existing indices onto the existing indices
      u <- unique(h)

      ## Transform the point set into the new indices
      Pt  <- P[u,]

      ## Transform the point correspondance mapping to the new index space  
      ht <- c()
      for (i in 1:length(h)) {
        ht[i] <- which(u == h[i])
      }

      ## DOESN'T WORK
      ## Form the inverse mapping from the existing indices to the new
      ## set of consecutive indices
      ## uinv <- c()
      ## uinv[u] <- 1:length(u)
      ## ht <- uinv[h[u]]

      ## Transform the triangulation to the new index space
      Tt  <- matrix(ht[T], ncol=3)

      ## Tansform the forward pointer into the new indices
      gft <- ht[gf]

      ## Determine H, the mapping from edges onto corresponding edges
      Cut <- matrix(ht[Cu], ncol=2)
      Cut <- t(apply(Cut, 1, sort))
      M <- nrow(Cut)
      H <- rep(0, M)
      for (i in 1:M) {
        if (!H[i]) {
          H[i] <- i
          for (j in i:M) {
            if (identical(Cut[i,], Cut[j,])) {
              H[j] <- i
            }
          }
        }
      }

      ## Form the mapping from a new set of consecutive edge indices
      ## onto the existing edge indices
      U <- unique(H)

      ## Transform the edge set into the new indices
      Cut <- Cut[U,]

      ## Transform the edge correspondance mapping to the new index space  
      Ht <- c()
      for (i in 1:length(H)) {
        Ht[i] <- which(U == H[i])
      }

      ## Create the lengths of the merged edges by averaging
      Lt <- c()
      for (k in 1:length(U)) {
        is <- which(Ht == k)
        ## if (length(is)>1) {
        ##   print(L[is])
        ## }
        Lt[k] <- mean(L[is])
      }

      ## Transform the rim set
      Rset <- order.Rset(self$ol$getRimSet(), self$ol$gf, self$ol$h)
      Rsett <- unique(ht[Rset])
      i0t <- ht[self$ol$i0]

      ## Create the symmetric connection set
      Ct <- rbind(Cut, Cut[,2:1])

      ## Matrix to map line segments onto the points they link
      ## Bt <- Matrix(0, nrow(Pt), nrow(Ct), sparse=TRUE)
      Bt <- matrix(0, nrow(Pt), nrow(Ct))
      for (i in 1:nrow(Ct)) {
        Bt[Ct[i,1],i] <- 1
      }

      self$Pt = Pt
      self$Tt = Tt
      self$Ct = Ct
      self$Cut = Cut
      self$Bt = Bt
      self$Lt = Lt
      self$ht = ht
      self$u = u
      self$U = U
      self$Rsett = Rsett
      self$i0t = i0t
      self$H = H
      self$Ht= Ht
    },
    ## This takes the mesh points from the flat outline and maps them to
    ## the curtailed sphere. It uses the area of the flat outline and
    ## \code{phi0} to determine the radius \code{R} of the sphere. It
    ## tries to get a good first approximation by using the function
    ## \code{\link{stretchMesh}}.
    ##
    ## @title Project mesh points in the flat outline onto a sphere
    ## @param r \code{Outline} object to which the following information
    ## has been added with \code{\link{mergePointsEdges}}:
    ## \describe{
    ## \item{\code{Pt}}{The mesh point coordinates.}
    ## \item{\code{Rsett}}{The set of points on the rim.}
    ## \item{\code{A.tot}}{The area of the flat outline.}}
    ## @return \code{reconstructedOutline} object containing the
    ## following extra information
    ## \item{\code{phi}}{Latitude of mesh points.}
    ## \item{\code{lmabda}}{Longitude of mesh points.}
    ## \item{\code{R}}{Radius of sphere.}
    ## @author David Sterratt
    ## @export
    projectToSphere = function() {
      Rsett <- self$Rsett
      i0t <- self$i0t
      A.tot <- self$ol$A.tot
      Cut <- self$Cut
      Lt <- self$Lt
      phi0 <- self$phi0
      lambda0 <- self$lambda0
      
      Nt <- nrow(self$Pt)
      Nphi <- Nt - length(Rsett)

      ## From this we can infer what the radius should be from the formula
      ## for the area of a sphere which is cut off at a latitude of phi0
      ## area = 2 * PI * R^2 * (sin(phi0)+1)
      R <- sqrt(A.tot/(2*pi*(sin(phi0)+1)))

      ## Find lengths between successive points on rim
      C <- matrix(NA, nrow(self$Pt), nrow(self$Pt))
      for (i in 1:nrow(Cut)) {
        C[Cut[i,1],Cut[i,2]] <- Lt[i]
        C[Cut[i,2],Cut[i,1]] <- Lt[i]
      }
      L.Rsett <- rep(NA, length(Rsett))
      for (i in 1:length(Rsett)) {
        L.Rsett[i] <- C[Rsett[i],Rsett[mod1(i+1, length(Rsett))]]
      }
      ## Check that this length matches the length computed from the AnnotatedOutline
      ## FIXME - this doesn't work for one retina - need to check why
      ## if (sum(L.Rsett) != sum(getRimLengths(r))) {
      ##  stop("Internal error: Mismatch in rim lengths")
      ## }
      ## Stretch mesh points to circle
      Ps <- stretchMesh(Cut, Lt, Rsett, circle(L=L.Rsett))
      x <- Ps[,1]
      y <- Ps[,2]
      phi <- -pi/2 + sqrt(x^2 + y^2)*(phi0+pi/2)
      phi[Rsett] <- phi0
      lambda <- atan2(y, x)
      lambda <- lambda - lambda[i0t] + lambda0

      self$phi <- phi
      self$lambda <- lambda
      self$R <- R
      self$phi0 <- phi0
      self$lambda0 <- lambda0
      self$Ps <- Ps
    },
    ## This function returns information about how edges on the sphere
    ## have been deformed from their flat state.
    ##
    ## @title Return strains edges are under in spherical retina
    ## @param r A \code{\link{ReconstructedOutline}} object
    ## @return A list containing two data frames \code{flat} and \code{spherical}. 
    ## Each data frame contains for each edge in the flat or spherical meshes:
    ## \item{\code{L}}{Length of the edge in the flat outline }
    ## \item{\code{l}}{Length of the corresponding edge on the sphere}
    ## \item{\code{strain}}{The strain of each connection}
    ## \item{\code{logstrain}}{The logarithmic strain of each connection}
    ## @author David Sterratt
    getStrains = function() {
      ## Original lengths in flattened outline is a vector with
      ## M elements, the number of rows of Cu
      L <- self$ol$L
      ## New lengths in reconstructed object is a vector with Mt < M
      ## elements, the number of rows of Cut
      lt <- compute.lengths(self$phi, self$lambda, self$Cut, self$R)
      ## For each connection in the flattened object, we want the length of
      ## the corresponding connection in the reconstructed object
      ## The mapping Ht achieves this
      l <- lt[self$Ht]
      stretch <- l/L
      strain <- stretch - 1
      logstrain <- log(stretch)

      ## Compute quantities in spherical retina too
      Lt <- self$Lt
      stretcht <- lt/Lt
      straint <- stretcht - 1
      logstraint <- log(stretcht)

      return(list(flat=
                    data.frame(L=L,  l=l,
                               strain=strain,  logstrain=logstrain),
                  spherical=
                    data.frame(L=Lt, l=lt,
                               strain=straint, logstrain=logstraint)))
    },
    ## Optimise the mapping from the flat outline to the sphere
    ##
    ## @title Optimise mapping
    ## @param r reconstructedOutline object
    ## @param alpha Area penalty scaling coefficient
    ## @param x0 Area penalty cut-off coefficient
    ## @param nu Power to which to raise area
    ## @param method Method to pass to \code{optim}
    ## @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
    ## @param dev.flat Device handle for plotting flatplot updates to. If
    ## \code{NA} don't make any flat plots
    ## @param dev.polar Device handle for plotting polar plot updates
    ## to. If \code{NA} don't make any polar plots.
    ## @param control Control argument to pass to \code{optim}
    ## @return reconstructedOutline object
    ## @author David Sterratt
    ## @export
    optimiseMapping = function(alpha=4, x0=0.5, nu=1, optim.method="BFGS",
                                plot.3d=FALSE, dev.flat=NA, dev.polar=NA,
                                control=list()) {
      phi <- self$phi
      lambda <- self$lambda
      R <- self$R
      phi0 <- self$phi0
      lambda0 <- self$lambda0
      Tt <- self$Tt
      A <- self$ol$A
      Cut <- self$Cut
      Ct <- self$Ct
      Lt <- self$Lt
      Bt <- self$Bt
      Rsett <- self$Rsett
      i0t <- self$i0t
      Nt <- nrow(self$Pt)
      Nphi <- Nt - length(Rsett)

      ## Optimisation and plotting
      opt <- list()
      opt$p <- c(phi[-Rsett], lambda[-i0t])
      opt$conv <- 1
      count <- 0
      while (opt$conv) {
        ## Optimise
        opt <- stats::optim(opt$p, E, gr=dE,
                            method=optim.method,
                            T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                            alpha=alpha,  N=Nt, x0=x0, nu=nu,
                            Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi,
                            verbose=FALSE, control=control)

        ## Report
        E.tot <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
                   alpha=alpha,  N=Nt, x0=x0, nu=nu,
                   Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
        E.l <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
                 alpha=0,  N=Nt, x0=x0, nu=nu,
                 Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)

        ft <- flipped.triangles(phi, lambda, Tt, R)
        nflip <- sum(ft$flipped)
        self$report(sprintf("E = %8.5f | E_L = %8.5f | E_A = %8.5f | %3d flippped triangles", E.tot, E.l, E.tot - E.l,  nflip))
        if (nflip & self$debug) {
          print(data.frame(rbind(id=which(ft$flipped),
                                 A=A[ft$flipped],
                                 a=ft$areas[ft$flipped])))
        }

        ## Decode p vector
        phi          <- rep(phi0, Nt)
        phi[-Rsett]  <- opt$p[1:Nphi]
        lambda       <- rep(lambda0, Nt)
        lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

        ## Plot
        if (plot.3d) {
          sphericalplot(list(phi=phi, lambda=lambda, R=R,
                             Tt=Tt, Rsett=Rsett, gb=self$ol$gb, ht=self$ol$ht),
                        datapoints=FALSE)
        }

        if (!is.na(dev.flat)) {
          dev.set(dev.flat)
          flatplot(self, grid=TRUE, strain=TRUE, mesh=FALSE, markup=FALSE,
                   datapoints=FALSE, landmarks=FALSE,
                   image=FALSE)
        }

        if (!is.na(dev.polar)) {
          ## Wipe any previous reconstruction of coordinates of pixels and feature sets
          self$ims <- NULL
          self$clearFeatureSets()
          dev.set(dev.polar)
          self$phi <- phi
          self$lambda <- lambda
          projection(self, mesh=TRUE, 
                     datapoints=FALSE, landmarks=FALSE,
                     image=FALSE)
        }
      }

      self$phi <- phi
      self$lambda <- lambda
      self$opt <- opt
      self$nflip <- sum(ft$flipped)
      self$E.tot <- E.tot
      self$E.l <- E.l
      self$mean.strain    <- mean(abs(self$getStrains()$spherical$strain))
      self$mean.logstrain <- mean(abs(self$getStrains()$spherical$logstrain))
    },
    ## Optimise the mapping from the flat outline to the sphere
    ##
    ## @title Optimise mapping
    ## @param r reconstructedOutline object
    ## @param alpha Area penalty scaling coefficient
    ## @param x0 Area penalty cut-off coefficient
    ## @param nu Power to which to raise area
    ## @param method Method to pass to \code{optim}
    ## @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
    ## @param dev.flat Device handle for plotting grid to
    ## @param dev.polar Device handle for plotting polar plot to
    ## @param ... Extra arguments to pass to \code{\link{fire}}
    ## @return reconstructedOutline object
    ## @author David Sterratt
    ## @export
    optimiseMappingCart  = function(alpha=4, x0=0.5, nu=1, method="BFGS",
                                 plot.3d=FALSE, dev.flat=NA, dev.polar=NA, ...) {
      phi <- self$phi
      lambda <- self$lambda
      R <- self$R
      phi0 <- self$phi0
      lambda0 <- self$lambda0
      Tt <- self$Tt
      A <- self$ol$A
      Cut <- self$Cut
      Ct <- self$Ct
      Lt <- self$Lt
      Bt <- self$Bt
      Rsett <- self$Rsett
      i0t <- self$i0t
      Nt <- nrow(self$Pt)
      Nphi <- Nt - length(Rsett)

      ## Optimisation and plotting
      opt <- list()
      opt$x <- sphere.spherical.to.sphere.cart(phi, lambda, R)
      opt$conv <- 1

      ## Compute "mass" for each node
      minL <- rep(Inf, nrow(self$Pt))
      for (i in 1:nrow(Cut)) {
        minL[Cut[i,1]] <- min(minL[Cut[i,1]], Lt[i])
        minL[Cut[i,2]] <- min(minL[Cut[i,2]], Lt[i])
      }
      m <- 1/minL
      m <- m/mean(m)
      count <- 50

      while (opt$conv && count) {
        ## Optimise
        opt <- fire(opt$x,
                    force=function(x) {Fcart(x, Ct, Lt, Tt, A, R, alpha, x0, nu)},
                    restraint=function(x) {Rcart(x, R, Rsett, i0t, phi0, lambda0)},
                    dt=1,
                    nstep=200,
                    m=m, verbose=TRUE, report=self$report, ...)
        count <- count - 1
        ## Report
        E.tot <- Ecart(opt$x, Cu=Cut, L=Lt, R=R, T=Tt, A=A,
                       alpha=alpha, x0=x0, nu=nu)
        E.l <- Ecart(opt$x, Cu=Cut, L=Lt, R=R, T=Tt, A=A,
                     alpha=0, x0=x0, nu=0)

        s <- sphere.cart.to.sphere.spherical(opt$x, R)
        phi <-    s[,"phi"]
        lambda <- s[,"lambda"]
        ft <- flipped.triangles(phi, lambda, Tt, R)
        nflip <- sum(ft$flipped)
        self$report(sprintf("E = %8.5f | E_L = %8.5f | E_A = %8.5f | %3d flippped triangles", E.tot, E.l, E.tot - E.l,  nflip))
        if (nflip) {
          print(data.frame(rbind(id=which(ft$flipped),
                                 A=A[ft$flipped],
                                 a=ft$areas[ft$flipped])))
        }

        ## Plot
        if (plot.3d) {
          sphericalplot(list(phi=phi, lambda=lambda, R=R,
                             Tt=Tt, Rsett=Rsett, gb=self$ol$gb, ht=self$ol$ht),
                        datapoints=FALSE)
        }

        if (!is.na(dev.flat)) {
          dev.set(dev.flat)
          flatplot(self, grid=TRUE, strain=TRUE, mesh=FALSE, markup=FALSE,
                   datapoints=FALSE, landmarks=FALSE,
                   image=FALSE)
        }

        if (!is.na(dev.polar)) {
          ## Wipe any previous reconstruction of coordinates of pixels and feature sets
          self$ims <- NULL
          self$clearFeatureSets()
          dev.set(dev.polar)
          self$phi <- phi
          self$lambda <- lambda
          projection(self, mesh=TRUE, 
                     datapoints=FALSE, landmarks=FALSE,
                     image=FALSE)
        }
      }

      self$phi <- phi
      self$lambda <- lambda
      self$opt <- opt
      self$nflip <- sum(ft$flipped)
      self$E.tot <- E.tot
      self$E.l <- E.l
      self$mean.strain    <- mean(abs(self$getStrains()$spherical$strain))
      self$mean.logstrain <- mean(abs(self$getStrains()$spherical$logstrain))
    },
    ## Transform an image into the reconstructed space. The four corner
    ## coordinates of each pixel are transformed into spherical
    ## coordinates and a mask matrix with the same dimensions as
    ## \code{im} is created. This has \code{TRUE} for pixels that should
    ## be displayed and \code{FALSE} for ones that should not.
    ##
    ## @title Transform an image into the reconstructed space
    ## @param r \code{reconstructedOutline} object
    ## @return \code{reconstructedOutline} object with extra elements
    ## \item{\code{ims}}{Coordinates of corners of pixes in spherical coordinates}
    ## \item{\code{immask}}{Mask matrix with same dimensions as image \code{im}}
    ## @author David Sterratt
    transformImage = function() {
      if (!is.null(self$ol$im)) {
        ## Need to find the *boundaries* of pixels
        N <- ncol(self$ol$im)
        M <- nrow(self$ol$im)

        ## Create grid coords of corners of pixels.  These run from the
        ## top left of the image down each column of the image.
        xs <- 0:N
        ys <- M:0
        ## x-coords of pixel corners, arranged in (N+1) by (M+1) grid
        ## Ditto for y-coords
        I <- cbind(as.vector(outer(ys*0, xs, FUN="+")),
                   as.vector(outer(ys, xs*0, FUN="+")))

        ## Find Barycentric coordinates of corners of pixels
        Ib <- tsearch(self$ol$getPoints()[,"X"], self$ol$getPoints()[,"Y"],
                      self$ol$T, I[,1], I[,2], bary=TRUE)
        rm(I)
        gc()
        ## Create mask depending on whether corners are in outline
        idx <- matrix(Ib$idx, M+1, N+1)
        self$immask <- (!is.na(idx[1:M    , 1:N    ]) &
                        !is.na(idx[1:M    , 2:(N+1)]) &
                        !is.na(idx[2:(M+1), 1:N    ]) &
                        !is.na(idx[2:(M+1), 2:(N+1)]))
        rm(idx)
        gc()
        ## Find 3D coordinates of mesh points
        Pc <- sph2cart(theta=self$lambda, phi=self$phi, r=1)
        self$ims <- bary2sph(Ib, self$Tt, Pc)
      }
    },
    ## Get coordinates of corners of pixels of image in spherical
    ## coordinates
    ## @param r \code{\link{ReconstructedOutline}} object
    ## @return Coordinates of corners of pixels in spherical coordinates
    ## @author David Sterratt
    ## @method getIms reconstructedOutline
    ## @export
    getIms = function(r) {
      if (is.null(self$ims)) {
        self$report("Transforming image...")
        ## Force garbage collection; not great practice, but this
        ## procedure is imemory intensive for large images
        gc()
        self$transformImage()
        gc()
      }
      return(self$ims)
    },
    ## @export
    getTearCoords = function(r) {
      Tss <- list()
      for (TF in self$ol$TFset) {
        ## Convert indices to the spherical frame of reference
        j <- self$ht[TF]
        Tss <- c(Tss, list(cbind(phi=self$phi[j], lambda=self$lambda[j])))
      }
      return(Tss)
    },
    getFeatureSet = function(type) {
      type <- paste0("Reconstructed", type)
      fs <- super$getFeatureSet(type)
      if (is.null(fs)) {
        self$reconstructFeatureSets()
        fs <- super$getFeatureSet(type)
      }
      return(fs)
    },
    reconstructFeatureSets = function() {
      self$featureSets <- lapply(self$ol$getFeatureSets(), function(x) x$reconstruct(self))
    }
  )
)
                           


##' Plot \code{\link{ReconstructedOutline}} object. This adds a mesh
##' of gridlines from the spherical retina (described by points
##' \code{phi}, \code{lambda} and triangulation \code{Tt} and cut-off
##' point \code{phi0}) onto a flattened retina (described by points
##' \code{P} and triangulation \code{T}).
##'
##' @title Flat plot of reconstructed outline
##' @param x \code{\link{ReconstructedOutline}} object
##' @param axt whether to plot axes
##' @param xlim x-limits
##' @param ylim y-limits
##' @param grid Whether or not to show the grid lines of
##' latitude and longitude
##' @param strain Whether or not to show the strain
##' @param ... Other plotting parameters
##' @method flatplot ReconstructedOutline
##' @author David Sterratt
##' @importFrom grDevices rainbow palette
##' @export
flatplot.ReconstructedOutline <- function(x, axt="n",
                                          xlim=NULL, ylim=NULL,
                                          grid=TRUE,
                                          strain=FALSE,
                                          ...) {
  ## NextMethod()
  flatplot(x$ol, ...)


  if (strain) {
    o <- x$getStrains()
    palette(rainbow(100))
    P <- x$ol$getPoints()
    scols <- strain.colours(o$flat$logstrain)
    Cu <- x$ol$Cu
    segments(P[Cu[,1],1], P[Cu[,1],2],
             P[Cu[,2],1], P[Cu[,2],2], col=round(scols))
  }
  
  ## Plot a gridline from the spherical retina (described by points phi,
  ## lambda and triangulation Tt) onto a flattened retina (described by
  ## points P and triangulation T). The gridline is described by a
  ## normal n to a plane and a distance to the plane. The intersection of
  ## the plane and the spehere is the gridline.
  get.gridline.flat <- function(P, T, phi, lambda, Tt, n, d, ...) {
    mu <- compute.intersections.sphere(phi, lambda, Tt, n, d)

    ## Take out rows that are not intersections. If a plane intersects
    ## one edge of a triangle and the opposing vertex, in the row
    ## corresponding to the triangle, there will be a 0, a 1 and a
    ## value between 0 and 1. We get rid of the 1 in the
    ## following. Triangles in which one line is in the plane have mu
    ## values 0, 1 and NaN; we want to include these.
    tri.int <- (rowSums((mu >= 0) & (mu <= 1), na.rm=TRUE) == 2) 
    ## | apply(mu, 1, function(x) setequal(x, c(0, 1, NaN))))
    if (any(tri.int)) {
      T  <- T[tri.int,,drop=FALSE]
      mu <- mu[tri.int,,drop=FALSE]

      ## Create a logical matrix of which points are involved in lines
      ## that interscect the plane.
      line.int <- (mu >= 0) & (mu < 1)
      ## If any element of mu contained a NaN, due to a line being in
      ## the plane, this should be set to false as the point opposite
      ## the NaN is not in the plane
      line.int[is.na(line.int)] <- FALSE

      ## Order rows so that the false indicator is in the third column
      T[!line.int[,2] ,] <- T[!line.int[,2], c(3,1,2)]
      mu[!line.int[,2],] <- mu[!line.int[,2],c(3,1,2)]
      T[!line.int[,1] ,] <- T[!line.int[,1], c(2,3,1)]
      mu[!line.int[,1],] <- mu[!line.int[,1],c(2,3,1)]

      P <- cbind(mu[,1] * P[T[,3],] + (1-mu[,1]) * P[T[,2],], 
                 mu[,2] * P[T[,1],] + (1-mu[,2]) * P[T[,3],])
                                        # suppressWarnings(segments(P1[,1], P1[,2], P2[,1], P2[,2], ...))
    } else {
      P <- matrix(0, nrow=1, ncol=4)
    }
    colnames(P) <- c("X1", "Y1", "X2", "Y2")
    return(P)
  }

  if (grid) {
    grid.int.minor <- 15
    grid.int.major <- 45
    grid.maj.col <- getOption("grid.maj.col")
    grid.min.col <- getOption("grid.min.col")

    phi0d <- x$phi0 * 180/pi

    P <- matrix(0, nrow=0, ncol=4)
    cols <- NULL
    Phis <- setdiff(seq(-90, phi0d, by=grid.int.minor), phi0d)
    Lambdas <- seq(0, 180-grid.int.minor, by=grid.int.minor)
    for (Phi in Phis) {
      if ((!(Phi %% grid.int.major) || Phi == phi0d)) {
        col <- grid.maj.col
      } else {
        col <- grid.min.col
      }
      P1 <- get.gridline.flat(x$ol$getPoints(), x$ol$T, x$phi, x$lambda, x$Tt,
                              c(0,0,1), sin(Phi*pi/180))
      cols <- c(cols, rep(col, nrow(P1)))
      P <- rbind(P, P1)
    }
    for (Lambda in Lambdas) {
      if (!(Lambda %% grid.int.major)) {
        col <- grid.maj.col
      } else {
        col <- grid.min.col
      }
      Lambda <- Lambda * pi/180
      P1 <- get.gridline.flat(x$ol$getPoints(), x$ol$T, x$phi, x$lambda, x$Tt,
                              c(sin(Lambda),cos(Lambda),0), 0)
      cols <- c(cols, rep(col, nrow(P1)))
      P <- rbind(P, P1)
    }
    if (nrow(P) > 0) {
      segments(P[,"X1"], P[,"Y1"], P[,"X2"], P[,"Y2"], col=cols)
    }
  }
}

##' Draw a projection of a \code{\link{ReconstructedOutline}}. This method sets up
##' the grid lines and the angular labels and draws the image.
##'
##' @title  Projection of a reconstructed outline
##' @param r \code{ReconstructedOutline} object
##' @param transform Transform function to apply to spherical coordinates
##' before rotation
##' @param axisdir Direction of axis (North pole) of sphere in
##' external space as matrix with column names \code{phi} (elevation)
##' and \code{lambda} (longitude).
##' @param projection Projection in which to display object,
##' e.g. \code{\link{azimuthal.equalarea}} or \code{\link{sinusoidal}}
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda} (longitude).
##' @param lambdalim Limits of longitude (in degrees) to display
##' @param philim Limits of latitude (in degrees) to display
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees
##' @param mesh If \code{TRUE}, plot mesh
##' @param grid Whether or not to show the grid lines of
##' latitude and longitude
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param colatitude If \code{TRUE} have radial labels plotted with
##' respect to colatitude rather than latitude
##' @param pole If \code{TRUE} indicate the pole with a "*"
##' @param image If \code{TRUE}, show the image
##' @param markup If \code{TRUE}, plot markup, i.e. reconstructed tears
##' @param add If \code{TRUE}, don't draw axes; add to existing plot.
##' @param max.proj.dim Maximum width of the image created in pixels
##' @param ... Graphical parameters to pass to plotting functions
##' @method projection ReconstructedOutline
##' @export
projection.ReconstructedOutline <- function(r,
                                            transform=identity.transform,
                                            axisdir=cbind(phi=90, lambda=0),
                                            projection=azimuthal.equalarea,
                                            proj.centre=cbind(phi=0, lambda=0),
                                            lambdalim=c(-180, 180),
                                            philim=c(-90, 90),
                                            labels=c(0, 90, 180, 270),
                                            mesh=FALSE,
                                            grid=TRUE,
                                            grid.bg="transparent",
                                            grid.int.minor=15,
                                            grid.int.major=45,
                                            colatitude=TRUE,
                                            pole=FALSE,
                                            image=TRUE,
                                            markup=TRUE,
                                            add=FALSE,
                                            max.proj.dim=getOption("max.proj.dim"),
                                            ...) {
  Call <- match.call(expand.dots=TRUE)
  plot.image <- image
  ## Compute grid lines

  ## Lines of latitude (parallels)

  ## Determine the major and minor parallels
  phis.maj <- seq(-90, 90, by=grid.int.major)
  phis.maj <- c(philim[1],
                phis.maj[(phis.maj > philim[1]) & (phis.maj < philim[2])],
                philim[2])
  phis.min <- seq(-90, 90, by=grid.int.minor)
  phis.min <- c(philim[0],
                phis.min[(phis.min > philim[1]) & (phis.min < philim[2])],
                philim[2])
  phis.min <- setdiff(phis.min, phis.maj)

  ## Longitudes at which to draw lines; the smaller the by interval,
  ## the smoother
  lambdas <- seq(lambdalim[1], lambdalim[2], by=1)

  ## Compute the minor and and major parallels to draw
  paras.min <- projection(pi/180*cbind(phi   =as.vector(outer(c(lambdas, NA)*0, phis.min, FUN="+")),
                                       lambda=as.vector(outer(c(lambdas, NA), phis.min*0, FUN="+"))),
                          proj.centre=pi/180*proj.centre)
  paras.maj <- projection(pi/180*cbind(phi   =as.vector(outer(c(lambdas, NA)*0, phis.maj, FUN="+")),
                                       lambda=as.vector(outer(c(lambdas, NA), phis.maj*0, FUN="+"))),
                          proj.centre=pi/180*proj.centre)

  ## Lines of longitude (meridians)

  ## Determine the major and minor parallels
  lambdas.maj <- c(rev(seq(0         , lambdalim[1], by=-grid.int.major)),
                   seq(grid.int.major, lambdalim[2], by= grid.int.major))
  lambdas.min <- c(rev(seq(0         , lambdalim[1], by=-grid.int.minor)),
                   seq(grid.int.minor, lambdalim[2], by= grid.int.minor))
  lambdas.min <- setdiff(lambdas.min, lambdas.maj)

  ## Lattidues at which to draw lines; the smaller the by interval,
  ## the smoother
  phis <- seq(philim[1], philim[2], by=1)

  ## Compute the minor and and major meridians to draw
  merids.min <- projection(pi/180*cbind(phi   =as.vector(outer(c(phis, NA), lambdas.min*0, FUN="+")),
                                        lambda=as.vector(outer(c(phis*0, NA), lambdas.min, FUN="+"))),
                           proj.centre=pi/180*proj.centre)
  merids.maj <- projection(pi/180*cbind(phi   =as.vector(outer(c(phis, NA), lambdas.maj*0, FUN="+")),
                                        lambda=as.vector(outer(c(phis*0, NA), lambdas.maj, FUN="+"))),
                           proj.centre=pi/180*proj.centre)

  ## Set up the plot region
  xlim <- range(na.omit(rbind(paras.min, paras.maj))[,"x"])
  ylim <- range(na.omit(rbind(paras.min, paras.maj))[,"y"])

  if (!add) {
    plot(NA, NA, xlim=xlim, ylim=ylim,# xaxs="i", yaxs="i",
         type = "n", axes = FALSE, xlab = "", ylab = "", asp=1,
         main=Call$main, bg=Call$bg)
  }

  ## Plot an image

  ## Get the spherical coordinates of the corners of pixels

  if (plot.image) {
    ims <- r$getIms()
    gc()
    if (!is.null(ims)) {
      ## Reconstitute image from stored values of phi and lambda
      ## coordinates of corners of pixels

      ## Get the size of the image
      M <- nrow(r$ol$im)
      N <- ncol(r$ol$im)

      ## Downsample the image by first selecting rows and columns to
      ## look at
      by <- ceiling(max(N, M)/max.proj.dim) # Number of pixels to merge
      Ms <- seq(1, M - (M %% by), by=by)
      Ns <- seq(1, N - (N %% by), by=by)

      ## Downsample the image
      ## This should not create a new version of the image
      im <- r$ol$im
      if (by > 1) {
        r$report("Downsampling image by factor of ", by)
        im <- r$ol$im[Ms, Ns]
      }

      ## Now need to do the more complex job of downsampling the matrix
      ## containing the coordinates of the corners of pixels
      if (by > 1) {
        r$report("Downsampling pixel corner spherical coordinates by factor of ", by)
        imsmask <- matrix(FALSE, M+1, N+1)
        imsmask[c(Ms, (max(Ms) + by)), c(Ns, (max(Ns) + by))] <- TRUE
        ims <- ims[imsmask,]
      }

      ## Convenience variables for the new image sizes
      M <- nrow(im)
      N <- ncol(im)

      ## Target size of of block, i.e. the number of pixels to
      ## transform in one go
      S <- 200*200 # 200*8000
      ## Number of blocks
      B <- ceiling(N*M/S)
      ## Numnber of columns in the first B-1 blocks
      C <- floor(N/B)
      ## In the Bth block there will be N-(B-1)*C columns
      ## print(paste("M =", M, "; N =", N, "; S =", S, "; B =", B, "; C =", C))
      
      ## Number of columns in block k
      Ck <- C
      for (k in 1:B) {
        r$report("Projecting block ", k, "/", B)
        ## Actual number of columns, since the last block may have a
        ## different number of chunks to C
        if (k == B) {
          Ck <- N - (B - 1)*C
        }
        ## Index of leftmost column of image matrix needed
        j0 <- (k - 1)*C + 1
        ## Index of rightmost column of image matrix needed
        j1 <- (k - 1)*C + Ck
      
        ## Transform the pixel coordinates and compute x and y positions
        ## of corners of pixels.
        ## inds <- (k - 1)*(M + 1)*C + (1:((M + 1)*(Ci + 1)))
        l0 <-  (j0 - 1)*(M + 1) + 1
        l1 <- j1*(M + 1) + (M + 1)
        ## print(paste("k =", k, "; Ck =", Ck, "; j0 =", j0, "; j1 =", j1,
        ##             "; l0 =", l0, "; l1 =", l1))
        inds <- l0:l1
        rc <- projection(
          rotate.axis(
            transform(ims[inds,],
                      phi0=r$phi0),
            axisdir*pi/180),
          lambdalim=lambdalim*pi/180,
          proj.centre=pi/180*proj.centre)

        xpos <- matrix(rc[,"x"], M + 1, Ck + 1)
        ypos <- matrix(rc[,"y"], M + 1, Ck + 1)

        ## Convert these to format suitable for polygon()
        impx <- rbind(as.vector(xpos[1:M    , 1:Ck    ]),
                      as.vector(xpos[1:M    , 2:(Ck+1)]),
                      as.vector(xpos[2:(M+1), 2:(Ck+1)]),
                      as.vector(xpos[2:(M+1), 1:Ck    ]),
                      NA)
        impy <- rbind(as.vector(ypos[1:M    , 1:Ck    ]),
                      as.vector(ypos[1:M    , 2:(Ck+1)]),
                      as.vector(ypos[2:(M+1), 2:(Ck+1)]),
                      as.vector(ypos[2:(M+1), 1:Ck   ]),
                      NA)
        ## Pixels outside the image should be masked. The mask is a matrix
        ## the same size as the image, containing TRUE for pixels that
        ## should be displayed and FALSE for those that should be
        ## masked. It is computed by finding the corners of the poly-pixel
        ## lie outwith the outline. These corners will have the coordinate
        ## NA.  print("sum(!is.na(colSums(impx[1:4,])))")
        ## print(sum(!is.na(colSums(impx[1:4,]))))
        immask <- matrix(!is.na(colSums(impx[1:4,])), M, Ck)

        ## We want to get rid of any poly-pixels that cross either end of
        ## the longitude range in a pseudocylindrical projection. A simple
        ## way of doing this is to say that if a pixel is very large,
        ## don't plot it.

        ## This code is chronically slow, hence replacing with the
        ## pmin/pmax version
        ## bigpx <- which(apply(impx[1:4,], 2,
        ##                      function(x) {max(x) - min(x)}) > 0.1*abs(diff(xlim)) |
        ##                apply(impy[1:4,], 2,
        ##                      function(y) {max(y) - min(y)}) > 0.1*abs(diff(ylim)))
        bigpx <- which((pmax(impx[1,], impx[2,], impx[3,], impx[4,]) -
                        pmin(impx[1,], impx[2,], impx[3,], impx[4,]) >
                        0.1*abs(diff(xlim))) |
                       (pmax(impy[1,], impy[2,], impy[3,], impy[4,]) -
                        pmin(impy[1,], impy[2,], impy[3,], impy[4,]) >
                        0.1*abs(diff(xlim))))
        immask[bigpx] <- FALSE

        ## Plot the polygon, masking as we go
        graphics::polygon(impx[,immask], impy[,immask],
                          col=im[1:M,(k-1)*C+(1:Ck)][immask],
                          border=im[1:M,(k-1)*C+(1:Ck)][immask])
      }
    }
  }

  ## Plot the mesh
  if (mesh) {
    Pt <- projection(rotate.axis(transform(cbind(phi=r$phi, lambda=r$lambda),
                                           phi0=r$phi0),
                                 axisdir*pi/180),
                     lambdalim=lambdalim*pi/180,
                     proj.centre=pi/180*proj.centre)
    trimesh(r$Tt, Pt, col="gray", add=TRUE)
  }
  
  grid.maj.col <- getOption("grid.maj.col")
  grid.min.col <- getOption("grid.min.col")

  ## Plot the grid
  if (grid) {
    ## Minor paralells and meridians
    lines(paras.min,  col=grid.min.col)
    lines(merids.min, col=grid.min.col)

    ## Major lines of latitude on top of all minor lines
    lines(paras.maj,  col=grid.maj.col)
    lines(merids.maj, col=grid.maj.col)

    ## Boundary of projection
    boundary <- projection("boundary")
    graphics::polygon(boundary[,"x"], boundary[,"y"], border="black")
  }

  ## Plot rim in visutopic space
  rs <- cbind(phi=r$phi0, lambda=seq(0, 2*pi, len=360))
  rs.rot <- rotate.axis(transform(rs, phi0=r$phi0), axisdir*pi/180)
  ## "Home" position for a cyclops looking ahead
  ## r$axisdir = cbind(phi=0, lambda=0)

  lines(projection(rs.rot, lambdalim=lambdalim*pi/180, lines=TRUE,
                   proj.centre=pi/180*proj.centre),
        col=getOption("TF.col"))

  ## Projection of pole
  if (pole) {
    oa.rot <- rotate.axis(transform(cbind(phi=-pi/2, lambda=0), phi0=r$phi0),
                          axisdir*pi/180)
    points(projection(oa.rot, proj.centre=pi/180*proj.centre),
           pch="*", col=getOption("TF.col"), cex=2)
  }

  ## Plot tears
  if (markup) {
    Tss <- r$getTearCoords()
    for (Ts in Tss) {
      ## Plot
      lines(projection(rotate.axis(transform(Ts, phi0=r$phi0),
                                   axisdir*pi/180),
                       lines=TRUE,
                       lambdalim=lambdalim*pi/180,
                       proj.centre=pi/180*proj.centre),
            col=getOption("TF.col"), lwd=Call$lwd, lty=Call$lty)
    }
  }

  if (grid) {  
    ## Longitude labels around rim - not on actual frame of reference!
    if (!is.null(labels)) {
      ## Longitudes (meridians) at which to plot at
      angles <- seq(0, by=2*pi/length(labels), len=length(labels))

      ## This is a nasty hack: we want to find out how far to plot the
      ## labels from the rim. We can't do this simply in terms of
      ## angles, so we have to find the right angular distance to give
      ## the desired fraction of the axes at which to plot the
      ## labels. This is done by this optimisation function.
      label.fax <- 0.02                   # Fraction of axis length from axes to plot labels
      opt <- stats::optimise(function(a) {
        rs0 <- cbind(phi=r$phi0,     lambda=angles[1])
        rs  <- cbind(phi=r$phi0 + a, lambda=angles[1])
        rc0 <- projection(rotate.axis(transform(rs0, phi0=r$phi0),
                                      axisdir*pi/180),
                          proj.centre=pi/180*proj.centre)
        rc  <- projection(rotate.axis(transform(rs, phi0=r$phi0),
                                      axisdir*pi/180),
                          proj.centre=pi/180*proj.centre)
        return((vecnorm(rc - rc0) - label.fax*abs(diff(xlim)))^2)
      }
     ,interval=c(1, 20)*pi/180)
      lambda.label.off <- opt$minimum

      ## Now plot the labels themselves. Phew!!
      rs <- cbind(phi=r$phi0 + lambda.label.off, lambda=angles)
      rc <- projection(rotate.axis(transform(rs, phi0=r$phi0), axisdir*pi/180),
                       proj.centre=pi/180*proj.centre)
      text(rc[,"x"], rc[,"y"], labels, xpd=TRUE)
    }

    ## Latitude Labels
    ## rlabels <- c(seq(philim[1], philim[2], by=grid.int.major))
    rlabels <- phis.maj
    rs <- cbind(phi=rlabels*pi/180, lambda=proj.centre[1,"lambda"])
    rc <- projection(rs, proj.centre=pi/180*proj.centre)
    text(rc[,"x"], rc[,"y"], rlabels + ifelse(colatitude, 90, 0),
         xpd=TRUE, adj=c(1, 1), col=grid.maj.col)
  }
}

##' @export
##' @importFrom graphics abline
lvsLplot.ReconstructedOutline <- function(r, ...) {
  Call <- match.call(expand.dots=TRUE)
  o <- r$getStrains()$spherical
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  cols <- strain.colours(o$logstrain)
  L <- o$L
  l <- o$l
  units <- r$ol$units
  xlab <- expression(paste("Flat ", italic(L)))
  ylab <- expression(paste("Reconstructed ", italic(l)))
  if (!is.na(units)) {
    xlab <- eval(substitute(expression(paste("Flat ", italic(L), " (", units1, ")")),
                           list(units1=units)))
    ylab <- eval(substitute(expression(paste("Reconstructed ", italic(l), " (",
                                        units1, ")")),
                            list(units1=units)))
  }
  
  
  plot(L, l, col=cols, pch=20,
       xlim=c(0, max(L, l)), ylim=c(0, max(L, l)),
       xlab=xlab, ylab=ylab,
       asp=1, main=Call$main)
  par(xpd=FALSE)
  abline(0, 1)
  abline(0, 0.75, col="blue")
  abline(0, 1.25, col="red")
  text(0.2*max(L), 0.2*max(L)*0.5, "25% compressed", col="blue",
               pos=4)
  text(0.75*max(L), 0.75*max(L)*1.25, "25% expanded", col="red",
               pos=2)
}

##' Draw a spherical plot of reconstructed outline. This method just
##' draws the mesh.
##'
##' @title Spherical plot of reconstructed outline
##' @param r \code{\link{ReconstructedOutline}} object
##' @param strain If \code{TRUE}, plot the strain
##' @param surf If \code{TRUE}, plot the surface
##' @param ... Other graphics parameters -- not used at present
##' @method sphericalplot reconstructedOutline
##' @author David Sterratt
##' @import rgl
##' @export
sphericalplot.reconstructedOutline <- function(r,
                                               strain=FALSE,
                                               surf=TRUE, ...) {
  NextMethod()
  
  ## FIXME: This needs to be looked at with a view to replacing
  ## functions in plots.R
  with(r, {
    ## Obtain Cartesian coordinates of points
    P <- sphere.spherical.to.sphere.cart(phi, lambda, R)

    if (surf) {
      ## Outer triangles
      fac <- 1.005
      triangles3d(matrix(fac*P[t(Tt[,c(2,1,3)]),1], nrow=3),
                  matrix(fac*P[t(Tt[,c(2,1,3)]),2], nrow=3),
                  matrix(fac*P[t(Tt[,c(2,1,3)]),3], nrow=3),
                  color="darkgrey", alpha=1)
      
      ## Inner triangles
      triangles3d(matrix(P[t(Tt),1], nrow=3),
                  matrix(P[t(Tt),2], nrow=3),
                  matrix(P[t(Tt),3], nrow=3),
                  color="white", alpha=1)
    }
    
    ## Plot any flipped triangles
    ft <- flipped.triangles(phi, lambda, Tt, R)
    with(ft, points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3],
                      col="blue", size=5))

    ## Shrink so that they appear inside the hemisphere
    fac <- 0.997
    rgl.lines(fac*rbind(P[ht[gb[gb]],1], P[ht[gb],1]),
              fac*rbind(P[ht[gb[gb]],2], P[ht[gb],2]),
              fac*rbind(P[ht[gb[gb]],3], P[ht[gb],3]),
              lwd=3, color=getOption("TF.col"))
    
    fac <- 1.006
    rgl.lines(fac*rbind(P[ht[gb[gb]],1], P[ht[gb],1]),
              fac*rbind(P[ht[gb[gb]],2], P[ht[gb],2]),
              fac*rbind(P[ht[gb[gb]],3], P[ht[gb],3]),
              lwd=3, color=getOption("TF.col"))

    if (strain) {
      o <- getStrains(r)
      palette(rainbow(100))
      scols <- strain.colours(o$spherical$logstrain)

      fac <- 0.999
      P1 <- fac*P[Cut[,1],]
      P2 <- fac*P[Cut[,2],]
      
      width <- 40
      ## Compute displacement vector to make sure that strips are
      ## parallel to surface of sphere
      d <- extprod3d(P1, P2-P1)
      d <- width/2*d/vecnorm(d)
      PA <- P1 - d
      PB <- P1 + d
      PC <- P2 + d
      PD <- P2 - d

      ## This is a ridiculously inefficient way of drawing the strain,
      ## but if you try presenting a color vector, it makes each line
      ## multi-coloured. It has taking HOURS of fiddling round to
      ## discover this! GRRRRRRRRRRRRRRRRR!
      for (i in 1:nrow(PA)) {
        quads3d(rbind(PA[i,1], PB[i,1], PC[i,1], PD[i,1]),
                rbind(PA[i,2], PB[i,2], PC[i,2], PD[i,2]),
                rbind(PA[i,3], PB[i,3], PC[i,3], PD[i,3]),
                color=round(scols[i]), alpha=1)
      }
    }
  })
}
