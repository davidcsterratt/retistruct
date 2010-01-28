source("M634-4.R")

par(mar=c(0,0,0,0))
plot.stitch(s, lwd=2)
with(t, trimesh(T, P, col="grey", lwd=2, add=TRUE))

dev.print(pdf, "../figures/M634-4-triangulated-stitched2.pdf", width=4, height=4)

with(s, plot.outline(P, gb, lwd=2))
with(as.list(c(t, m, p)), plot.gridlines.flat(P, T, phi, lambda, Tt, phi0, lwd=2))

dev.print(pdf, "../figures/M634-4-initial-proj2.pdf", width=4, height=4)

with(as.list(c(t, m, p)), plot.retina(phi, lambda, R, Tt, Rsett))
view3d(0, -20, zoom=0.7)
rgl.postscript("../figures/M634-4-initial-proj-3d2.pdf", fmt="pdf")

with(s, plot.outline(P, gb, lwd=2))
with(as.list(c(t, m, p)), plot.gridlines.flat(P, T, r$phi, r$lambda, Tt, phi0, lwd=2))
dev.print(pdf, "../figures/M634-4-final-proj2.pdf", width=4, height=4)

with(as.list(c(t, m, p)), plot.retina(r$phi, r$lambda, R, Tt, Rsett))
view3d(0, -20, zoom=0.7)
rgl.postscript("../figures/M634-4-final-proj-3d2.pdf", fmt="pdf")
