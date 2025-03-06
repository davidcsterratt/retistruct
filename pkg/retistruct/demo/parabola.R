dataset <- file.path(system.file(package = "retistruct"), "extdata", "smi32-csv")
## Read in unmorphed retina
a <- retistruct.read.markup(retistruct.read.dataset(dataset))

## Read in morphed retina and plot
focal.length <- 50                      # Focal length of parabola
ap <- retistruct.read.markup(retistruct.read.dataset(morph.dataset.to.parabola(orig.dataset=dataset, f=focal.length)))
## tp <- TriangulatedDataset(TriangulatedOutline(ap))
depthplot3D(ap)

## Reconstruct both
rp <- retistruct.reconstruct(ap)
r  <- retistruct.reconstruct(a)

message("Residual energy of original retina:", r$E.l)
message("Residual energy of transformed retina:", rp$E.l)

oldpar <- par(no.readonly=TRUE) # Save graphics parameters before plotting
par(mfcol=c(2, 3))
par(mar=c(1, 1, 1, 1))
flatplot(r, mesh=TRUE, main="Original retina")
flatplot(rp, mesh=TRUE, main="Parabolically transformed retina")
projection(r, projection=azimuthal.conformal, datapoint.contours=FALSE, main="Conformal projection of original retina", mesh=TRUE)
projection(rp, projection=azimuthal.conformal, datapoint.contours=FALSE, main="Conformal projection of transformed retina", mesh=TRUE)
par(mar=c(5, 4, 1, 1), mgp=c(1.5, 0.5, 0))
lvsLplot(r, main="Strain plot of original retina")
lvsLplot(rp, main="Strain plot of transformed retina")

par(oldpar) # Restore graphics parameters
