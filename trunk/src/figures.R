source("M634-4.R")

## r <- optimise.mapping(p, m, t, s, E0.A=100)
r <- optimise.mapping(p, m, t, s, E0.A=exp(3), k.A=1)
p1 <- p
p1$phi <- r$phi
p1$lambda <- r$lambda
r <- optimise.mapping(p1, m, t, s, E0.A=exp(10), k.A=20)

par(mar=c(0,0,0,0))
plot.stitch(s, lwd=2)
points(P.red[,1], P.red[,2], col="red", pch=".", cex=5)
points(P.green[,1], P.green[,2], col="green", pch=".", cex=5)
with(t, trimesh(T, P, col="grey", lwd=1, add=TRUE))
with(s, plot.outline(P, gb, lwd=2, col="orange", add=TRUE))

dev.print(pdf, "../figures/M634-4-triangulated-stitched2.pdf", width=4, height=4)

with(s, plot.outline(P, gb, lwd=2, col="orange"))
points(P.red[,1], P.red[,2], col="red", pch=".", cex=5)
points(P.green[,1], P.green[,2], col="green", pch=".", cex=5)
with(as.list(c(t, m, p)), plot.gridlines.flat(P, T, phi, lambda, Tt, phi0, lwd=1))
with(s, plot.outline(P, gb, lwd=2, col="orange", add=TRUE))
dev.print(pdf, "../figures/M634-4-initial-proj2.pdf", width=4, height=4)

with(as.list(c(t, m, p)), plot.retina(phi, lambda, R, Tt, Rsett))
plot.outline.retina(p$phi, p$lambda, p$R*1.01, s$gb, m$ht, lwd=5, color="orange")
with(as.list(c(t, m, p)), plot.cell.bodies(phi, lambda, R, Tt, cb.red, col="red"))
with(as.list(c(t, m, p)), plot.cell.bodies(phi, lambda, R, Tt, cb.green, col="green"))
view3d(0, -20, zoom=0.7)
rgl.postscript("../figures/M634-4-initial-proj-3d2.pdf", fmt="pdf")

with(s, plot.outline(P, gb, lwd=2))
points(P.red[,1], P.red[,2], col="red", pch=".", cex=5)
points(P.green[,1], P.green[,2], col="green", pch=".", cex=5)
with(as.list(c(t, m, p)), plot.gridlines.flat(P, T, r$phi, r$lambda, Tt, phi0, lwd=1))
with(s, plot.outline(P, gb, lwd=2, col="orange", add=TRUE))
dev.print(pdf, "../figures/M634-4-final-proj2.pdf", width=4, height=4)

with(as.list(c(t, m, p)), plot.retina(r$phi, r$lambda, R, Tt, Rsett))
plot.outline.retina(r$phi, r$lambda, p$R*1.01, s$gb, m$ht, lwd=5, color="orange")
with(as.list(c(t, m, p)), plot.cell.bodies(phi, lambda, R, Tt, cb.red, col="red"))
with(as.list(c(t, m, p)), plot.cell.bodies(phi, lambda, R, Tt, cb.green, col="green"))
view3d(0, -20, zoom=0.7)
rgl.postscript("../figures/M634-4-final-proj-3d2.pdf", fmt="pdf")


## Basic
##r = solve.mapping.momentum(p, m, t, s, E0.A=0, dt=1E-4, nstep=40000, Rexp=1, verbose=FALSE)

## Rexp 2 - doesn't work with this step size
## r1 = solve.mapping.momentum(p, m, t, s, E0.A=0, dt=1E-4, nstep=40000, Rexp=2, verbose=FALSE)





