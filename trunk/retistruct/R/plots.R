##
## Utilities for plotting functions
## 

## Generate colours for strain plots
strain.colours <- function(x) {
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  col <- x/log(0.75)*35 + 35
  col[col<1] <- 1
  col[col>70] <- 70
  return(col)
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
