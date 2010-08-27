##
## Plotting functions
## 

## Each function takes a list containing members from the orginal data
## or derived from the reconstruction procedure
##
## General format for function name is
##
## plot.<structure>.<view>
##
## where <structure> is one of
## - outline
## - stitch
## - sphere
## - gridlines
## - datapoints
## - strain
##
## and <view> is one of
## - flat [default]
## - spherical
## - polar
## - polararea
##

##
## Flat plots
##

## plot.outline.flat(P, gb)
##
## Plot outline of retina given set of outline points P and backwards
## pointer gb
plot.outline.flat <- function(P, gb, add=FALSE, axt="n", ...) {
  s <- which(!is.na(gb))                # source index
  d <- na.omit(gb)                      # destination index
  if (!add) {
    par(mar=c(1.4, 1.4, 1, 1), mgp=c(2, 0.2, 0), tcl=-0.2)
    plot(P[s,1], P[s,2], pch=".", xaxt=axt, yaxt=axt, xlab="", ylab="",
         bty="n")
  }
  segments(P[s,1], P[s,2], P[d,1], P[d,2], ...)
}

## plot.stitch.flat(P, s)
##
## Plot stitch given set of outline points stitch information s
plot.stitch.flat <- function(s, add=FALSE, ...) {
  with(s, {
    if (!add) plot.outline.flat(P, gb, ...)
    points(P[VF,], col="red", pch="+")
    points(P[VB,], col="orange", pch="+")
    points(P[V0,], col="cyan", pch="+")
    for (TF in TFset) {
      lines(P[TF,], col="red", ...)
    }
    for (TB in TBset) {
      lines(P[TB,], col="orange", ...)
    }
    for (j in 1:length(h)) {
      if (h[j] != j) {
        lines(P[c(j, h[j]),], col="blue", ...)
      }
    }
    
    for (j in 1:length(hf)) {
      if (hf[j] != j) {
        lines(P[c(j, hf[j]),], col="green", ...)
      }
    }
  })
}

## plot.gridline.flat(P, T, phi, lambda, Tt, n, d)
##
## Plot a gridline from the spherical retina (described by points phi,
## lambda and triangulation Tt) onto a flattened retina (described by
## points P and triangulation T). The gridline is described by a
## normal n to a plane and a distance to the plane. The intersection of
## the plane and the spehere is the gridline.
plot.gridline.flat <- function(P, T, phi, lambda, Tt, n, d, ...) {
  mu <- compute.intersections.sphere(phi, lambda, Tt, n, d)

  ## Take out rows that are not intersections
  tri.int <- (rowSums((mu >=0) & (mu <=1)) == 2)

  if (any(tri.int)) {
    T  <- T[tri.int,,drop=FALSE]
    mu <- mu[tri.int,,drop=FALSE]

    line.int <- (mu >=0) & (mu <=1)
    
    ## Order rows so that the false indicator is in the third column
    T[!line.int[,2] ,] <- T[!line.int[,2], c(3,1,2)]
    mu[!line.int[,2],] <- mu[!line.int[,2],c(3,1,2)]
    T[!line.int[,1] ,] <- T[!line.int[,1], c(2,3,1)]
    mu[!line.int[,1],] <- mu[!line.int[,1],c(2,3,1)]

    P1 <- mu[,1] * P[T[,3],] + (1-mu[,1]) * P[T[,2],]
    P2 <- mu[,2] * P[T[,1],] + (1-mu[,2]) * P[T[,3],]
    segments(P1[,1], P1[,2], P2[,1], P2[,2], ...)
  }
}

## plot.gridlines.flat(P, T, phi, lambda, Tt, phi0)
##
## Plot a mesh of gridlines from the spherical retina (described by
## points phi, lambda and triangulation Tt and cutoff point phi0) onto
## a flattened retina (described by points P and triangulation T).
plot.gridlines.flat <- function(P, T, phi, lambda, Tt, phi0,
                                Phis=(-8:9)*pi/18, Lambdas=(0:17)*pi/18,
                                grid.int.minor=15, grid.int.major=45,
                                grid.col="gray", ...) {
  Phis <- seq(-90, phi0, by=grid.int.minor)
  Lambdas <- seq(0, 180-grid.int.minor, by=grid.int.minor)
  for (Phi in Phis) {
    if ((!(Phi %% grid.int.major) || Phi == phi0)) {
      col <- "black"
    } else {
      col <- grid.col
    }
    plot.gridline.flat(P, T, phi, lambda, Tt, c(0,0,1), sin(Phi*pi/180), col=col, ...)
  }
  for (Lambda in Lambdas) {
    if (!(Lambda %% grid.int.major)) {
      col <- "black"
    } else {
      col <- grid.col
    }
    Lambda <- Lambda * pi/180
    plot.gridline.flat(P, T, phi, lambda, Tt, c(sin(Lambda),cos(Lambda),0), 0, col=col, ...)
  }
}

## Function to plot data points on flat outline
##
## Arguments
## Ds     - list of sets of datapoints, in which name of each set
##          is the colour in which to plot
##
plot.datapoints.flat <- function(Ds, ...) {
  for(col in names(Ds)) {
    points(Ds[[col]][,1], Ds[[col]][,2], col=col, pch=20,cex=0.5, ...)
  }
}

## Function to plot landmarks on the flat outline
##
## Arguments
## Ss     - list of segments
##
plot.landmarks.flat <- function(Ss, ...) {
  if (length(Ss) > 0) {
    for(i in 1:length(Ss)) {
      lines(Ss[[i]][,1], Ss[[i]][,2], ...)
    }
  }
}

## Generate colours for strain plots
strain.colours <- function(x) {
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  col <- x/log(0.75)*35 + 35
  col[col<1] <- 1
  col[col>70] <- 70
  return(col)
}

## Function to plot the fractional change in length of connections 
plot.strain.flat <- function(r) {
  o <- compute.strain(r)
  cols <- strain.colours(log(o$strain))
  with(r, 
       segments(P[Cu[,1],1], P[Cu[,1],2],
                P[Cu[,2],1], P[Cu[,2],2], col=cols))
}

## Function to plot the fractional change in length of connections 
plot.l.vs.L <- function(r) {
  o <- compute.strain(r)
  op <- par()["mar"]
  par(mar=c(4.5, 4.5, 0.5,0.5))
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  cols <- strain.colours(log(o$strain))
  with(o, plot(L, l, col=cols, pch='.', cex=5,
               xlim=c(0, max(L, l)), ylim=c(0, max(L, l)),
               xlab="Length on flattened object",
               ylab="Length on reconstructed object", asp=1))
  par(xpd=FALSE)
  abline(0, 1)
  abline(0, 0.75, col="blue")
  abline(0, 1.25, col="red")
  with(o, text(0.7*max(L), 0.7*max(L)*0.75, "25% compressed", col="blue",
               pos=4))
  with(o, text(0.75*max(L), 0.75*max(L)*1.25, "25% expanded", col="red",
               pos=2))
  par(op)
}

##
## Spherical plots
##

## Function to plot the mesh describing the reconstructed hemisphere
## in 3D
##
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Tt     - triagulation
## Rsett  - members of rim set
##
plot.sphere.spherical <- function(phi, lambda, R, Tt, Rsett) {
  ## Now plot this in 3D space....
  x <- R*cos(phi)*cos(lambda) 
  y <- R*cos(phi)*sin(lambda)
  z <- R*sin(phi)
  P <- cbind(x, y, z)
  rgl.clear()
  rgl.bg(color="white")

  ## Outer triangles
  triangles3d(matrix(1.01*x[t(Tt[,c(2,1,3)])], nrow=3),
              matrix(1.01*y[t(Tt[,c(2,1,3)])], nrow=3),
              matrix(1.01*z[t(Tt[,c(2,1,3)])], nrow=3),
              color="darkgrey", alpha=1)
  
  ## Inner triangles
  triangles3d(matrix(x[t(Tt)], nrow=3),
              matrix(y[t(Tt)], nrow=3),
              matrix(z[t(Tt)], nrow=3),
              color="white", alpha=1)

  ## Plot any flipped triangles
  ft <- flipped.triangles(phi, lambda, Tt, R)
  with(ft, points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3],
                    col="blue", size=5))
}

## Function to plot outline in 3D
## 
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## gb     - outline pointer
## h      - outline correspondences
##
plot.outline.spherical <- function(phi, lambda, R, gb, h, ...) {
  ## Obtain Cartesian coordinates of points
  Pc <- sphere.spherical.to.sphere.cart(phi, lambda, R)

  ## Shrink so that they appear inside the hemisphere
  P <- Pc*0.99
  rgl.lines(rbind(P[h[gb[gb]],1], P[h[gb],1]),
            rbind(P[h[gb[gb]],2], P[h[gb],2]),
            rbind(P[h[gb[gb]],3], P[h[gb],3]),
             ...)
  
  P <- Pc*1.001
  rgl.lines(rbind(P[h[gb[gb]],1], P[h[gb],1]),
            rbind(P[h[gb[gb]],2], P[h[gb],2]),
            rbind(P[h[gb[gb]],3], P[h[gb],3]),
             ...)
}

## Function to plot data points on a sphere
##
## It assumes that plot.sphere.spherical has been called already
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Tt     - triagulation
## cb     - object returned by tsearch containing information on the
##          triangle in which a cell body is found and its location
##          within that triangle in barycentric coordinates
## size   - size of the points to plot
## color  - colour of the points to plot
plot.datapoints.spherical <- function(phi, lambda, R, Tt, cb, size=R/10, color="red") {
  ## Obtain Cartesian coordinates of points
  cc <- datapoints.sphere.cart(phi, lambda, R, Tt, cb)
  
  ## Plot
  ## shade3d( translate3d( cube3d(col=color), cc[,1], cc[,2], cc[,3]))
  ## rgl.spheres(cc[,1], cc[,2], cc[,3], radius, color=color)
  ## points3d(cc[,1], cc[,2], cc[,3], size=size, color=color)
  ## cc <- cc * 1.01
  ## points3d(cc[,1], cc[,2], cc[,3], size=size, color=color)

  ## Custom code required to plot triangles
  ax1 <- 1/sqrt(apply(cc[,1:2]^2, 1, sum)) * cbind(-cc[,2], cc[,1], 0)
  ## print(ax1)
  
  ax2 <- extprod3d(cc, ax1)
  ax2 <- ax2/sqrt(apply(ax2^2, 1, sum))
  ##print(ax2)

  ##  print(dot(ax1, ax2))
  
  v1 <- cc + size *  ax1/2
  v2 <- cc + size * (-ax1/4 + sqrt(3)/4*ax2)
  v3 <- cc + size * (-ax1/4 - sqrt(3)/4*ax2)

  inmag <- 0.99
  outmag <- 1.02
  
  x <- rbind(v2[,1], v1[,1], v3[,1])
  y <- rbind(v2[,2], v1[,2], v3[,2])
  z <- rbind(v2[,3], v1[,3], v3[,3])
  triangles3d(inmag*x, inmag*y, inmag*z, color=color)

  x <- rbind(v1[,1], v2[,1], v3[,1])
  y <- rbind(v1[,2], v2[,2], v3[,2])
  z <- rbind(v1[,3], v2[,3], v3[,3])
  triangles3d(outmag*x, outmag*y, outmag*z, color=color)
}

##
## Polar plots
##

## Function to set up polar plot
## 
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## gb     - outline pointer
## h      - outline correspondences
##
plot.polar <- function(phi0=40,
                       show.grid=TRUE, grid.col="gray", grid.bg="transparent",
                       grid.int.minor=15, grid.int.major=45,
                       radial.labels.major=c("N", "", "D", "", "T", "", "V", "")) {
  par(mar=c(1, 1, 1, 1))
  grid.pos <- c(seq(-90, phi0, by=grid.int.minor), phi0)
  maxlength <- diff(range(grid.pos))
  plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
       type = "n", axes = FALSE, xlab = "", ylab = "", asp=1)

  ## Plot the grid
  ## Tangnential lines
  angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
  if (show.grid) {
    for (i in seq(length(grid.pos), 1, by = -1)) {
      xpos <- cos(angles) * (grid.pos[i] - grid.pos[1])
      ypos <- sin(angles) * (grid.pos[i] - grid.pos[1])
      if (((grid.pos[i] %% grid.int.major) == 0) || (i == length(grid.pos))) {
        col <- "black"
      } else {
        col <- grid.col
      }
      polygon(xpos, ypos, border = col, col = grid.bg)
    }
  }

  ## Radial lines
  angles <- seq(0, 180-grid.int.minor, by = grid.int.minor)
  col <- rep(grid.col, length(angles))
  col[(angles %% grid.int.major) == 0] <- "black"
  angles <- angles * pi/180
  xpos <- cos(angles) * maxlength
  ypos <- sin(angles) * maxlength
  segments(xpos, ypos, -xpos, -ypos, col=col)

  ## Tangential Labels
  labels <- c(seq(-90, phi0, by=grid.int.major), phi0)
  label.pos <- labels - min(grid.pos)
  text(label.pos, -maxlength/15, labels)

  ## Radial Labels
  angles <- seq(0, 360-grid.int.major, by=grid.int.major)*pi/180
  xpos <- cos(angles) * maxlength * 1.05
  ypos <- sin(angles) * maxlength * 1.05
  text(xpos, ypos, radial.labels.major)
}

## Function to plot the outline in polar coordinates
plot.outline.polar <- function(r, ...) {
  ## Look through forward tears
  for (TF in r$TFset) {
    ## Convert indicies to the spherical frame of reference
    j <- r$ht[TF]
    ## Plot
    x <- with(r, cos(lambda[j]) * ((phi[j] * 180/pi) + 90))
    y <- with(r, sin(lambda[j]) * ((phi[j] * 180/pi) + 90))
    lines(x, y, ...)
  }
}

## Function to plot cell bodies in spherical coordinates on a polar plot
##
## Arguments
## Dss     - list of sets of datapoints, each of which is defined by a matrix
##           with columns "phi" and "lambda".
##           The name of each element in the list is the colour in which the
##           datapoints in each set are plotted.
##
plot.datapoints.polar <- function(Dss, pch=".", ...) {
  for (i in 1:length(Dss)) {
    phis    <- Dss[[i]][,"phi"]
    lambdas <- Dss[[i]][,"lambda"]
    xpos <- cos(lambdas) * ((phis * 180/pi) + 90)
    ypos <- sin(lambdas) * ((phis * 180/pi) + 90)      
    points(xpos, ypos, col=names(Dss)[i], pch=pch, ...)
  }
}

## Function to plot landmarks on the polar plot
##
## Arguments
## Sss     - list of landmarks, each of which is defined by a matrix
##           with columns "phi" and "lambda"
##
plot.landmarks.polar <- function(Sss, ...) {
  if (length(Sss) > 0) {
    for (i in 1:length(Sss)) {
      phi    <- Sss[[i]][,"phi"]
      lambda <- Sss[[i]][,"lambda"]
      x <- cos(lambda) * ((phi * 180/pi) + 90)
      y <- sin(lambda) * ((phi * 180/pi) + 90)
      lines(x, y, ...)
    }
  }
}

## Function to plot segments in the polar plot
plot.segments.polar <- function(phi0, lambda0, phi1, lambda1, ...) {
  x0 <- cos(lambda0) * ((phi0 * 180/pi) + 90)
  y0 <- sin(lambda0) * ((phi0 * 180/pi) + 90)
  x1 <- cos(lambda1) * ((phi1 * 180/pi) + 90)
  y1 <- sin(lambda1) * ((phi1 * 180/pi) + 90)
  segments(x0, y0, x1, y1, ...)
}

## 
## Polar area plots
##

## Function to plot cell bodies in spherical coordinates on a polar plot
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Tt     - triagulation
## cbs    - list of objects returned by tsearch containing information on the
##          triangle in which a cell body is found and its location
##          within that triangle in barycentric coordinates
## phi0   - lattitude of the rim in radians
## cols   - colour of points to plot for each object in cbs
plot.datapoints.polararea <- function(phi, lambda, R, Tt, cbs, phi0, cols="red",
                                   pch=".", ...) {
  plot(NA, NA, xlim=c(-2,2), ylim=c(-2, 2))
  for (i in 1:length(cbs)) {
    cs <- datapoints.sphere.spherical(phi, lambda, R, Tt, cbs[[i]])
    ## Turn into polar coordinates, shifting round by 90 degress for plotting
    lambdas <- cs$lambda+pi/2
    p <- polar.to.cart(spherical.to.polar.area(cs$phi), lambdas)
    points(p[,"x"], p[,"y"], pch=pch, col=cols[i], ...)

    ## Compute mean and plot
    m <- sphere.mean.sphere(cs$phi, lambdas)
    p <- polar.to.cart(spherical.to.polar.area(m[1]), m[2])
    points(p[,"x"], p[,"y"], col=cols[i], pch="+", ...)
  }
  ## Draw circular grid
  dl <- 2*pi/90
  lambdas <- seq(dl, 2*pi, by=dl)
  phi.degs <- seq(-80, phi0*180/pi, by=10)
  rs <- spherical.to.polar.area(phi.degs*pi/180)
  polygon(rbind(outer(cos(lambdas), rs), NA),
          rbind(outer(sin(lambdas), rs), NA),  col=NA, border="grey")

  ## Draw axes and label
  axis(side=1, pos=0, at=c(-max(rs), 0, 1, max(rs)), labels=c(NA, -90, 0, phi0*180/pi))
  axis(side=2, pos=0, at=c(-max(rs), max(rs)), labels=c(NA, NA))
  text(2, 0, "N")
  text(-2, 0, "T")
  text(0, 2, "D")
  text(0, -2, "V")
  
}
