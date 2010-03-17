source("M634-4.R")

## Optimisation procedure
## r <- optimise.mapping(p, m, t, s, E0.A=100)
r <- optimise.mapping(p, m, t, s, E0.A=exp(3), k.A=1)
p1 <- p
p1$phi <- r$phi
p1$lambda <- r$lambda
r <- optimise.mapping(p1, m, t, s, E0.A=exp(10), k.A=20)

## Figure of triangulated and stitched retina
par(mar=c(0,0,0,0))
plot.stitch(s, lwd=2)
points(P.red[,1], P.red[,2], col="red", pch=".", cex=5)
points(P.green[,1], P.green[,2], col="green", pch=".", cex=5)
with(t, trimesh(T, P, col="grey", lwd=1, add=TRUE))
with(s, plot.outline(P, gb, lwd=2, col="orange", add=TRUE))

dev.print(pdf, "../figures/M634-4-triangulated-stitched2.pdf", width=4, height=4)

## Figure of initial projection in flat retina space
with(s, plot.outline(P, gb, lwd=2, col="orange"))
points(P.red[,1], P.red[,2], col="red", pch=".", cex=5)
points(P.green[,1], P.green[,2], col="green", pch=".", cex=5)
with(as.list(c(t, m, p)), plot.gridlines.flat(P, T, phi, lambda, Tt, phi0, lwd=1))
with(s, plot.outline(P, gb, lwd=2, col="orange", add=TRUE))
dev.print(pdf, "../figures/M634-4-initial-proj2.pdf", width=4, height=4)

## Figure of initial projection in folded retinal space
with(as.list(c(t, m, p)), plot.retina(phi, lambda, R, Tt, Rsett))
plot.outline.retina(p$phi, p$lambda, p$R*1.01, s$gb, m$ht, lwd=5, color="orange")
with(as.list(c(t, m, p)), plot.cell.bodies(phi, lambda, R, Tt, cb.red, col="red"))
with(as.list(c(t, m, p)), plot.cell.bodies(phi, lambda, R, Tt, cb.green, col="green"))
view3d(0, -20, zoom=0.6)
rgl.postscript("../figures/M634-4-initial-proj-3d2.pdf", fmt="pdf")
rgl.snapshot("../figures/M634-4-initial-proj-3d2.png", top=TRUE)

## Figure of final projection in flat retinal space
with(s, plot.outline(P, gb, lwd=2))
points(P.red[,1], P.red[,2], col="red", pch=".", cex=5)
points(P.green[,1], P.green[,2], col="green", pch=".", cex=5)
with(as.list(c(t, m, p)), plot.gridlines.flat(P, T, r$phi, r$lambda, Tt, phi0, lwd=1))
with(s, plot.outline(P, gb, lwd=2, col="orange", add=TRUE))
dev.print(pdf, "../figures/M634-4-final-proj2.pdf", width=4, height=4)

## Figure of final projection in folded retinal space
with(as.list(c(t, m, p)), plot.retina(r$phi, r$lambda, R, Tt, Rsett))
plot.outline.retina(r$phi, r$lambda, p$R*1.01, s$gb, m$ht, lwd=5, color="orange")
with(as.list(c(t, m, p)), plot.cell.bodies(r$phi, r$lambda, R, Tt, cb.red, col="red"))
with(as.list(c(t, m, p)), plot.cell.bodies(r$phi, r$lambda, R, Tt, cb.green, col="green"))
view3d(0, -20, zoom=0.6)
rgl.postscript("../figures/M634-4-final-proj-3d2.pdf", fmt="pdf")
rgl.snapshot("../figures/M634-4-final-proj-3d2.png", top=TRUE)

## Figure of final projection in polar plot
with(as.list(c(t, m, p)), plot.cell.bodies.polar(r$phi, r$lambda, R, Tt, list(cb.red, cb.green), phi0, cols=c("red", "green"), cex=5))
dev.print(pdf, "../figures/M634-4-final-proj-polar2.pdf", width=4, height=4)

## Basic
##r = solve.mapping.momentum(p, m, t, s, E0.A=0, dt=1E-4, nstep=40000, Rexp=1, verbose=FALSE)

## Rexp 2 - doesn't work with this step size
## r1 = solve.mapping.momentum(p, m, t, s, E0.A=0, dt=1E-4, nstep=40000, Rexp=2, verbose=FALSE)





