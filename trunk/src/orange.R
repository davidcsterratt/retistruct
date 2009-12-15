source("fold-sphere2.R")
library("pixmap")

im <- read.pnm("../data/20091214191056-small.ppm")
# plot(im)

## Read in data
P.raw <- read.table("../data/20091214191056-small.txt")
P <- cbind(P.raw[,1], 520-P.raw[,2])

tearmat <- rbind(c(5,   1, 11),
                 c(21, 16, 32),
                 c(29, 25, 31),
                 c(52, 39, 65))
colnames(tearmat) <- c("apex", "end1", "end2")


## Stitching
s <- stitch.retina(P, tearmat)
plot.stitch(s)

t <- make.triangulation(s, d=30)
with(t, trimesh(T, P, col="black"))

m <- merge.points(t, s)
## Plot stiched retina in 2D (messy)
## trimesh(Tt, Pt, col="black")

## Plotting
plot(P)
with(s, plot.outline(P, gb))

p <- project.to.sphere(m, t, phi0=45*pi/180)

## Initial plot in 3D space
plot.retina(p$phi, p$lambda, p$R, m$Tt, m$Rsett)

## Optimisation
f <- solve.mapping(p, m, t, s, dt=1E-5, nstep=1000)

plot(im)
plot.gridlines.flat(t$P, t$T, f$phi, f$lambda, m$Tt, p$phi0,
                                Phis=(-1:2)*pi/4, Lambdas=(0:1)*pi/2, lwd=2)
