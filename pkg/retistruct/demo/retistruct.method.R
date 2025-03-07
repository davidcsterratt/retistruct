library(rgl)
oldpar <- par(no.readonly=TRUE) # Save graphics parameters before plotting
## Set up a 3x3 grid for plotting
par(mfrow=c(3, 3))
par(mar=c(0.5, 0.5, 0.5, 0.5))

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM509/R-CONTRA")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Make this a left eye to help with orientation of points
o$side="Left"

## Plot of raw data. Axes are reversed to improve comparison with
## polar plot later
flatplot(o, markup=FALSE)
mtext("A", adj=0, font=2, line=-0.9)

## Plot the annotation
flatplot(o, datapoints=FALSE, landmarks=FALSE)
mtext("B", adj=0, font=2, line=-0.9)

## Set up fixed point
o$lambda0 <- 0

## In v0.6 and below, the code below could be used to triangulate and
## stitch the outline. In v0.7,
## RetinalReconstructedOutline$loadOutline() carries out these steps.

## Initial triangulation (with 500 points)
# n <- 500
# t <- TriangulatedOutline$new(o, n=n)

## Stitching
# s <- StitchedOutline(t)

## Triangulate again, to take into account points added by stitching
# m <- TriangulatedOutline(s, n=n,
#                          suppress.external.steiner=TRUE)

## Merge the points that have been stitched
# m <- mergePointsEdges(m)

## Make a rough projection to a sphere
# m <- projectToSphere(m)

r <- NULL
r <- RetinalReconstructedOutline$new()
r$loadOutline(o, debug=FALSE)

## Plot triangulation and stitching
flatplot(r, datapoints=FALSE, landmarks=FALSE, markup=FALSE)
mtext("C", adj=0, font=2, line=-0.9)

## Plot the initial gridlines in 2D
par(mfg=c(3, 3))
flatplot(r, grid=TRUE,
         datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
         stitch=FALSE, strain=TRUE)
mtext("Dii", adj=0, font=2, line=-0.9)

## Plot the intial projection in 3D
par(mfg=c(2, 3))
plot.new()
mtext("Di", adj=0, font=2, line=-0.9)

sphericalplot(r, strain=TRUE, datapoints=FALSE)
view3d(zoom=0.7)
## Save to SVG
# rgl.postscript("initial-projection.svg", "svg")

r$reconstruct(plot.3d=FALSE)

## Plot the final projection in 3D and on the grid
par(mfg=c(2, 2))
plot.new()
mtext("Ei", adj=0, font=2, line=-0.9)

sphericalplot(r, strain=TRUE, datapoints=FALSE)
## Save to SVG
# rgl.postscript("final-projection.svg", "svg")

par(mfg=c(3, 2))
flatplot(r, grid=TRUE,
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE)
mtext("Eii", adj=0, font=2, line=-0.9)

## Plot data in polar coordinates and flattend retina
par(mfg=c(2, 1))
projection(r, datapoints=TRUE, landmarks=TRUE, datapoint.contours=FALSE)
mtext("Fi", adj=0, font=2, line=-0.9)

par(mfg=c(3, 1))
flatplot(r, grid=TRUE,
          datapoints=TRUE, landmarks=TRUE, mesh=FALSE, markup=FALSE,
          stitch=FALSE)
mtext("Fii", adj=0, font=2, line=-0.9)

## Save to PDF
# dev.print(pdf, file="retistruct-method.pdf", width=6.83, height=6.83)

par(oldpar) # Restore graphics parameters
