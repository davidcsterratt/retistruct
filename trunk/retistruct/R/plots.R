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
  o <- getStrains(r)
  cols <- strain.colours(o$logstrain)
  with(r, 
       segments(P[Cu[,1],1], P[Cu[,1],2],
                P[Cu[,2],1], P[Cu[,2],2], col=cols))
}

## Function to plot the fractional change in length of connections 
plot.l.vs.L <- function(r) {
  o <- getStrains(r)
  op <- par()["mar"]
  par(mar=c(4.5, 4.5, 0.5,0.5))
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  cols <- strain.colours(o$logstrain)
  with(o, plot(L, l, col=cols, pch=20,
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

## Function to plot data points on a sphere
##
## It assumes that plot.sphere.spherical has been called already
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Dsc    - structure containing locations of datapoints in spherical
##          cartesian coordinates
## size   - size of the points to plot
plot.datapoints.spherical <- function(phi, lambda, R, Dsc, D.cols, size=R/10) {
  for(col in names(Dsc)) {
    Dc <- Dsc[[col]]                      # Cartesian coordinates of points

    ## Code required to plot triangles

    ## Find two axes that are orthogonal to line from the origin to
    ## the datapoint
    
    ## Find axis in z=0 plane that is orthogonal to projection of
    ## datapoint onto that plane
    ax1 <- 1/sqrt(apply(Dc[,1:2]^2, 1, sum)) * cbind(-Dc[,2], Dc[,1], 0)

    ## Find axis that is orthogonal to the plane of axis 1 and the
    ## datapoint
    ax2 <- extprod3d(Dc, ax1)
    ax2 <- ax2/sqrt(apply(ax2^2, 1, sum))

    ## Create the verticies of an equillateral triangle to plot
    v1 <- Dc + size *  ax1/2
    v2 <- Dc + size * (-ax1/4 + sqrt(3)/4*ax2)
    v3 <- Dc + size * (-ax1/4 - sqrt(3)/4*ax2)

    ## Plot the triangle inside and outside the sphere
    inmag <- 0.99
    outmag <- 1.02
  
    x <- rbind(v2[,1], v1[,1], v3[,1])
    y <- rbind(v2[,2], v1[,2], v3[,2])
    z <- rbind(v2[,3], v1[,3], v3[,3])
    triangles3d(inmag*x, inmag*y, inmag*z, color=D.cols[[col]])

    x <- rbind(v1[,1], v2[,1], v3[,1])
    y <- rbind(v1[,2], v2[,2], v3[,2])
    z <- rbind(v1[,3], v2[,3], v3[,3])
    triangles3d(outmag*x, outmag*y, outmag*z, color=D.cols[[col]])
  }
}

##
## Polar plots
##

## Put text on the polar plot
text.polar <- function(text) {
  mtext(text, 1, adj=1, line=-0.1)
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
plot.datapoints.polararea <- function(phi, lambda, R, Tt, cbs, phi0,
                                      D.cols, cols="red",
                                      pch=".", ...) {
  plot(NA, NA, xlim=c(-2,2), ylim=c(-2, 2))
  for (i in 1:length(cbs)) {
    cs <- datapoints.sphere.spherical(phi, lambda, R, Tt, cbs[[i]])
    ## Turn into polar coordinates, shifting round by 90 degress for plotting
    lambdas <- cs$lambda+pi/2
    p <- polar.to.cart(spherical.to.polar.area(cs$phi), lambdas)
    points(p[,"x"], p[,"y"], pch=pch, col=D.cols[[cols[i]]], ...)

    ## Compute mean and plot
    m <- sphere.mean.sphere(cs$phi, lambdas)
    p <- polar.to.cart(spherical.to.polar.area(m[1]), m[2])
    points(p[,"x"], p[,"y"], col=D.cols[[cols[i]]], pch="+", ...)
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
