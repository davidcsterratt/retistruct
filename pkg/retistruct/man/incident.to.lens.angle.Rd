\name{incident.to.lens.angle}
\alias{incident.to.lens.angle}
\title{Convert incident angle to angle on lens}
\usage{
  incident.to.lens.angle(psi, anp, S)
}
\arguments{
  \item{psi}{Incident angles}

  \item{anp}{Anterior nodal point on x-axis}

  \item{S}{Lens system}
}
\value{
  Angles on lens
}
\description{
  Convert incident angle to angle on lens
}
\examples{
S <- S.mouse.ss.23
draw.system(S, xoff=-2)
## imfile <- "eye23.png"
## im <- as.raster(readPNG(imfile))
scale <- 0.0086
xoff <- -0.52
yoff <- 0.31
## rasterImage(im, xoff, yoff-scale*nrow(im)/2, xoff+scale*ncol(im), yoff+scale*nrow(im)/2)
draw.system(S, add=TRUE, col="yellow")
anp <- 1
psi <- 0.1
psi <- 0.2
## Trace ray angles using nodal point at 0.2 - not too bad an approximation
## trace.ray.angle(0.9, 0.2, S)
## All that's needed is to get the angle at which the ray hits the final lens.
## Then we have a mapping from psi to phi.
## r <- trace.ray(c(-500, 101.5, -0.2), S)
angles <- seq(0, 90, by=10)
lens.angles <- 180/pi*incident.to.lens.angle(angles*pi/180, 0.2, S)
}
\author{
  David Sterratt
}

